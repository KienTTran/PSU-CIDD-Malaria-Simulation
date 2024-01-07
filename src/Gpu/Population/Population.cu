//
// Created by kient on 12/31/2023.
//

#include <curand_kernel.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/sort.h>
#include "Population.cuh"
#include "Population/Population.h"
#include "Model.h"
#include "Core/Config/Config.h"
#include "MDC/ModelDataCollector.h"
#include "Spatial/SpatialModel.hxx"
#include "Core/Random.h"
#include "Gpu/Utils.cuh"
#include "Gpu/Random.cuh"
#include "Population/Properties/PersonIndexByLocationMovingLevel.h"

GPU::Population::Population() {
}

void GPU::Population::init() {
}
/*
 * Calculate number of circulations from each location
 * this if from poisson distribution with mean = popsize_residence_by_location * circulation_percent
 * d_n_circulations_from_locations is number of circulations from each location
 * */
__global__ void calculate_circulation_number(int n_locations,
                                             int* d_popsize_residence_by_location,
                                             double circulation_percent,
                                             curandState* d_state,
                                             int* d_n_circulations_from_locations) {
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    curandState local_state = d_state[thread_index];
    for(int index = thread_index; index < n_locations; index += stride) {
        if (d_popsize_residence_by_location[index] == 0 || circulation_percent == 0.0) return;
        auto poisson_means = d_popsize_residence_by_location[index] * circulation_percent;
        d_n_circulations_from_locations[index] = curand_poisson(&d_state[index], poisson_means);
    }
    d_state[thread_index] = local_state;
}

/*
 * Calculate relative out movement to destination locations
 * this is get_v_relative_out_movement_to_destination in CPU v4.0 version
 * Output is relative out movement from each location to each destination location
 * So 9 locations will have 9*9 = 81 values
 * */
__global__ void calculate_circulation_probabilities(int n_locations,
                                                    int* d_popsize_residence_by_location,
                                                    double* d_spatial_model_parameters,
                                                    double* d_spatial_model_travels,
                                                    int* d_district_vector,
                                                    double* d_distance_vector,
                                                    double* d_relative_outmovement,
                                                    int* d_from_indices,
                                                    int* d_target_indices) {
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations*n_locations; index += stride) {
        int from_location = index / n_locations;
        int target_location = index % n_locations;
        if(d_distance_vector[from_location * n_locations + target_location] == 0.0) return;
        double distance = d_distance_vector[from_location * n_locations + target_location];
        double kernel = pow(1 + (distance / d_spatial_model_parameters[2]), (-d_spatial_model_parameters[1]));
        double probability = pow(d_popsize_residence_by_location[from_location], d_spatial_model_parameters[0]) * kernel;
        probability = probability / (1 + d_spatial_model_travels[from_location] + d_spatial_model_travels[target_location] );
        if (d_district_vector[from_location] == static_cast<int>(d_spatial_model_parameters[3]) &&
            d_district_vector[target_location] == static_cast<int>(d_spatial_model_parameters[3])) {
            probability /= d_spatial_model_parameters[4];
        }
        d_relative_outmovement[index] = probability;
        d_from_indices[index] = from_location;
        d_target_indices[index] = target_location;
    }
}

void GPU::Population::calculate_circulate_locations(int n_locations,ThrustTVectorDevice<int> &d_n_circulations_from_locations,ThrustTVectorDevice<double> &d_relative_outmovement_from_target_,
                                         ThrustTVectorDevice<int> &d_all_location_from_indices,ThrustTVectorDevice<int> &d_all_location_target_indices) {
    //Has to get pointer to device otherwise it will copy vector from host to device
    d_ce_popsize_residence_by_location = Model::DATA_COLLECTOR->popsize_residence_by_location_gpu();
    d_ce_spatial_model_parameters = Model::CONFIG->spatial_model()->getSpatialModelParameters();
    d_ce_spatial_model_travels = Model::CONFIG->spatial_model()->getSpatialModelTravels();
    d_ce_spatial_districts = Model::CONFIG->h_spatial_districts;
    d_ce_spatial_distances = Model::CONFIG->h_spatial_distances;

    //All probabilities because thrust run all arrays at the same time. If use 1 array then values are overwritten.
    d_n_circulations_from_locations.resize(n_locations,0);

    //Get circulations by location
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int block_size = (d_n_circulations_from_locations.size() + n_threads - 1)/n_threads;
    calculate_circulation_number<<<block_size,n_threads>>>(n_locations,
                                                           thrust::raw_pointer_cast(d_ce_popsize_residence_by_location.data()),
                                                           Model::CONFIG->circulation_info().circulation_percent,
                                                           Model::GPU_RANDOM->d_states,
                                                           thrust::raw_pointer_cast(d_n_circulations_from_locations.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());

    //Get outmovement probabilities
    d_relative_outmovement_from_target_.resize(n_locations*n_locations,0.0);
    d_all_location_from_indices.resize(n_locations*n_locations,0);
    d_all_location_target_indices.resize(n_locations*n_locations,0);
    n_threads = Model::CONFIG->gpu_config().n_threads;
    block_size = (d_relative_outmovement_from_target_.size() + n_threads - 1)/n_threads;
    calculate_circulation_probabilities<<<block_size,n_threads>>>(Model::CONFIG->number_of_locations(),
                                                                  thrust::raw_pointer_cast(d_ce_popsize_residence_by_location.data()),
                                                                  thrust::raw_pointer_cast(d_ce_spatial_model_parameters.data()),
                                                                  thrust::raw_pointer_cast(d_ce_spatial_model_travels.data()),
                                                                  thrust::raw_pointer_cast(d_ce_spatial_districts.data()),
                                                                  thrust::raw_pointer_cast(d_ce_spatial_distances.data()),
                                                                  thrust::raw_pointer_cast(d_relative_outmovement_from_target_.data()),
                                                                  thrust::raw_pointer_cast(d_all_location_from_indices.data()),
                                                                  thrust::raw_pointer_cast(d_all_location_target_indices.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());
}

/*
 * Calculate moving level density at each location, output size is n_location*n_moving_level
 * This value is for each destination locations
 * */
__global__ void calculate_moving_level_density_kernel(int n_locations,
                                               int n_moving_levels,
                                               thrust::tuple<int,int> *d_circulation_indices,
                                               int* d_popsize_by_moving_level,
                                               double* d_moving_level_value,
                                               double* d_moving_level_density) {
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations*n_moving_levels; index += stride) {
        int circulate_loc_index = index / n_moving_levels;
        int from_location = thrust::get<0>(d_circulation_indices[circulate_loc_index]);
        int moving_level = index % n_moving_levels;
        if(d_popsize_by_moving_level[from_location*n_moving_levels+moving_level] == 0 || d_moving_level_value[moving_level] == 0.0) return;
        d_moving_level_density[index] = d_popsize_by_moving_level[from_location*n_moving_levels+moving_level] * d_moving_level_value[moving_level];
    }
}

void GPU::Population::calculate_moving_level_density(ThrustT2TupleVectorDevice<int,int> d_circulation_indices,ThrustTVectorDevice<double> &d_moving_level_density) {
    if(d_circulation_indices.size() == 0){
        return;
    }
    d_moving_level_density.resize(d_circulation_indices.size()*Model::CONFIG->circulation_info().number_of_moving_levels);
    d_ce_popsize_by_moving_level = Model::CONFIG->h_popsize_by_moving_level;
    d_ce_moving_level_value = Model::CONFIG->circulation_info().v_moving_level_value;
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int block_size = (d_moving_level_density.size() + n_threads - 1)/n_threads;
    calculate_moving_level_density_kernel<<<block_size,n_threads>>>(d_circulation_indices.size(),
                                                             Model::CONFIG->circulation_info().number_of_moving_levels,
                                                             thrust::raw_pointer_cast(d_circulation_indices.data()),
                                                             thrust::raw_pointer_cast(d_ce_popsize_by_moving_level.data()),
                                                             thrust::raw_pointer_cast(d_ce_moving_level_value.data()),
                                                             thrust::raw_pointer_cast(d_moving_level_density.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());
}

template <typename T>
struct notZero : public thrust::unary_function<T,bool> {
    __host__ __device__
    bool operator()(T x)
    {
        return x != 0;
    }
};

struct copyNotZero : public thrust::unary_function<unsigned int,bool>{
    __host__ __device__
    bool operator()(unsigned int x) {
        return x != 0;
    }
};

template <typename T>
struct isOne : public thrust::unary_function<T,bool> {
    __host__ __device__
    bool operator()(T x)
    {
        return x == 1;
    }
};

struct circulateLess{
    __host__ __device__
    bool operator()(const thrust::tuple<int,int,int,unsigned int>& t1, thrust::tuple<int,int,int,unsigned int>& t2)
    {
        if(t1.get<3>() < t2.get<3>())
            return true;
        if(t1.get<3>() > t2.get<3>())
            return false;
        return t1.get<3>() < t2.get<3>();
    }
};

/*
 * Get 3 vectors of from_location, target_location and moving_level for removing zero from
 * d_circulation_indices and d_n_circulations_all_loc_ml
 * */
__global__ void extract_locations_and_moving_levels(int n_locations,
                                          int n_moving_levels,
                                          thrust::tuple<int,int> *d_circulation_indices,
                                          unsigned int *d_n_circulations_all_loc_ml,
                                          int* d_from_indices,
                                          int* d_target_indices,
                                          int* d_moving_levels){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations*n_moving_levels; index += stride) {
        if(d_n_circulations_all_loc_ml[index] == 0) return;
        int circulate_index = index / n_moving_levels;
        d_from_indices[index] = thrust::get<0>(d_circulation_indices[circulate_index]);
        d_target_indices[index] = thrust::get<1>(d_circulation_indices[circulate_index]);
        d_moving_levels[index] = index % n_moving_levels;
    }
}

/*
 * Combine 2 vectors to 1 vectors, for sorting
 * */
__global__ void zip_location_indices_and_n_circulations(int size,
                                          thrust::tuple<int,int,int>* d_circulations_all_loc_ml_indices_no_zero,
                                          unsigned int* d_n_circulations_all_loc_ml_no_zero,
                                          thrust::tuple<int,int,int,unsigned int>* d_circulate_all_loc_ml_today){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < size; index += stride) {
        if(d_n_circulations_all_loc_ml_no_zero[index] == 0) return;
        d_circulate_all_loc_ml_today[index] = thrust::make_tuple(thrust::get<0>(d_circulations_all_loc_ml_indices_no_zero[index]),
                                                                     thrust::get<1>(d_circulations_all_loc_ml_indices_no_zero[index]),
                                                                     thrust::get<2>(d_circulations_all_loc_ml_indices_no_zero[index]),
                                                                     d_n_circulations_all_loc_ml_no_zero[index]);
    }
}

/*
 * Get random person index at each moving level in each location in all locations
 * Note that this is just filling 1 person index for each moving level in each location
 * so n_persons data in d_circulate_all_loc_ml_today is not used but will be embedded in d_circulate_person_indices_today
 * */
__global__ void fill_circulate_person_indices(int n_locations,
                                             int n_moving_levels,
                                             curandState *d_state,
                                             thrust::tuple<int,int,int,unsigned int> *d_circulate_all_loc_ml_today,
                                             int * d_popsize_by_loc_ml,
                                             thrust::tuple<int,int,int,unsigned int,int> *d_circulate_person_indices_today){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    curandState local_state = d_state[thread_index];
    for(int index = thread_index; index < n_locations; index += stride) {
        if(thrust::get<3>(d_circulate_all_loc_ml_today[index]) == 0) return;
        int from_location = thrust::get<0>(d_circulate_all_loc_ml_today[index]);
        int target_location = thrust::get<1>(d_circulate_all_loc_ml_today[index]);
        int moving_level = thrust::get<2>(d_circulate_all_loc_ml_today[index]);
        unsigned int n_persons = thrust::get<3>(d_circulate_all_loc_ml_today[index]);
        int size = d_popsize_by_loc_ml[from_location*n_moving_levels+moving_level];
        if(size == 0 || n_persons == 0) return;
        /*
         * Random uniform to get person index
         * To get result same as gsl, using casting method, which is [0,n-1]
         * ceiling methos is [1,n]
         * https://github.com/nglee/so_answers/blob/master/cuda/170426/kernel.cu
         * */
//        printf("kernel %d from_location: %d, moving_level: %d, popsize_loc_ml: %d size: %d\n",index,
//               from_location,moving_level,popsize_loc_ml,d_n_circulations_all_loc_ml_no_zero[index]);
        if(n_persons == 1){
            int p_index = curand_uniform(&local_state) * size;
            d_circulate_person_indices_today[index] = thrust::make_tuple(from_location,target_location,moving_level,n_persons,p_index);
//        printf("kernel %d from_location: %d, moving_level: %d, popsize_loc_ml: %d, size: %d, p_index: %d\n",index,from_location,moving_level,
//               popsize_loc_ml,d_n_circulations_all_loc_ml_no_zero[index],p_index);
        }
        else{
            d_circulate_person_indices_today[index] = thrust::make_tuple(from_location,target_location,moving_level,n_persons,-1);
        }
    }
    d_state[thread_index] = local_state;
}

/*
 * To speed up circulation process using GPU, doing following steps:
 * 1. Calculate number of circulations and leavers, same as CPU
 * 2. Filter out location with zero circulation
 * 3. Parallel multinomial in all non-zero circulation locations,
 * each location select number of leavers in each moving level
 * 4. Sort the result of multinomial, so all location and moving level with 1 person will be filled first on GPU
 * 5. For the rest, fill on CPU
 * */
void GPU::Population::perform_circulation_event() {
    auto tp_start = std::chrono::high_resolution_clock::now();

    /*
     * Calculate probability of leaving location in all locations (n_location*n_location)
     * Also get indices from and to arrays
     * */
    Model::GPU_POPULATION->calculate_circulate_locations(Model::CONFIG->number_of_locations(),
                                                         d_ce_n_circulations_from_locations,
                                                         d_ce_relative_outmovement_from_target,
                                                         d_ce_all_location_from_indices, d_ce_all_location_target_indices);

    /*
     * Calculate number of leavers in all locations
     * */
    ThrustTVectorDevice<unsigned int> d_num_leavers_from_target_(d_ce_relative_outmovement_from_target.size(), 0);
    Model::GPU_RANDOM->random_multinomial(Model::CONFIG->number_of_locations(),
                                          Model::CONFIG->number_of_locations(),
                                          d_ce_n_circulations_from_locations,
                                          d_ce_relative_outmovement_from_target,
                                          d_num_leavers_from_target_);
    size_t no_zero_size = thrust::count_if(d_num_leavers_from_target_.begin(), d_num_leavers_from_target_.end(), notZero<unsigned int>());
    /*
     * Remove zero values in d_num_leavers_from_target_ and d_n_circulations_from_locations
     * d_circulations_indices_no_zero index is not location index, its index is index of n_locations*n_locations
     * thrust::get<0>(d_circulations_indices_no_zero[i]) is from location index
     * thrust::get<1>(d_circulations_indices_no_zero[i]) is to location index
     * */
    ThrustT2TupleVectorDevice<int,int> d_circulations_indices_no_zero(no_zero_size);

    /*
     * Remove zero values in d_num_leavers_from_target_ and d_n_circulations_from_locations
     * scan d_num_leavers_from_target_ and copy non-zero from & to locations to d_circulations_indices_no_zero
     * This is to reduce compute time
     */
    auto loc_index_begin = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.begin(), d_ce_all_location_target_indices.begin()));
    auto loc_index_end = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.end(), d_ce_all_location_target_indices.end()));
    thrust::copy_if(thrust::device,
                     loc_index_begin,
                     loc_index_end,
                     d_num_leavers_from_target_.begin(),
                     d_circulations_indices_no_zero.begin(),
                     copyNotZero());
    thrust::device_vector<unsigned int>::iterator nend = thrust::remove(thrust::device,d_num_leavers_from_target_.begin(),d_num_leavers_from_target_.end(),0);
    ThrustTVectorDevice<unsigned int> d_num_leavers_from_target_no_zero(no_zero_size);
    thrust::copy(d_num_leavers_from_target_.begin(),nend,d_num_leavers_from_target_no_zero.begin());

    /*
     * d_num_leavers_from_target_no_zero is non-zero leavers in all locations
     * */
    int total_leavers = thrust::reduce(thrust::device,d_num_leavers_from_target_no_zero.begin(),d_num_leavers_from_target_no_zero.end());
//    LOG_IF(total_leavers == 0, DEBUG) << "[Population] Update population circulation GPU total_leavers = 0";
    if(total_leavers == 0) return;

    /*
     * Calculate moving level density at each location, output size is n_location*n_moving_level
     * n_location is no-zero indices
     */
    ThrustTVectorDevice<double> d_moving_level_density;
    Model::GPU_POPULATION->calculate_moving_level_density(d_circulations_indices_no_zero,d_moving_level_density);

    ThrustTVectorDevice<unsigned int> d_n_circulations_all_loc_ml(d_moving_level_density.size(),0);
    Model::GPU_RANDOM->random_multinomial(d_circulations_indices_no_zero.size(),
                                         Model::CONFIG->circulation_info().number_of_moving_levels,
                                         d_num_leavers_from_target_no_zero,
                                         d_moving_level_density,
                                         d_n_circulations_all_loc_ml);
    /*
     * d_circulations_indices_no_zero and d_n_circulations_all_loc_ml not the same size
     */
//    ThrustTVectorHost<unsigned int> h_n_circulations_all_loc_ml = d_n_circulations_all_loc_ml;
//    ThrustT2TupleVectorHost<int,int> h_circulations_indices_no_zero = d_circulations_indices_no_zero;
//    for(int i = 0; i < h_n_circulations_all_loc_ml.size(); i++){
//        int circulate_loc_index = i / Model::CONFIG->circulation_info().number_of_moving_levels;
//        int from_location = thrust::get<0>(h_circulations_indices_no_zero[circulate_loc_index]);
//        int target_location = thrust::get<1>(h_circulations_indices_no_zero[circulate_loc_index]);
//        int moving_level = i % Model::CONFIG->circulation_info().number_of_moving_levels;
//        if(h_n_circulations_all_loc_ml[i] == 0) continue;
//        printf("%d from %d to %d moving level %d size %d\n",i,from_location,target_location,moving_level,
//               h_n_circulations_all_loc_ml[i]);
//    }
//    printf("\n");

    /*
     * Remove zero values in d_n_circulations_all_loc_ml
     * First, extract from and to locations from d_n_circulations_all_loc_ml to 2 vectors
     */
    d_ce_all_location_from_indices.resize(d_n_circulations_all_loc_ml.size());
    d_ce_all_location_target_indices.resize(d_n_circulations_all_loc_ml.size());
    d_ce_all_moving_levels.resize(d_n_circulations_all_loc_ml.size());
    if(d_n_circulations_all_loc_ml.size() == 0) return;
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int block_size = (d_n_circulations_all_loc_ml.size() + n_threads - 1)/n_threads;
    extract_locations_and_moving_levels<<<block_size,n_threads>>>(d_circulations_indices_no_zero.size(),
                                                      Model::CONFIG->circulation_info().number_of_moving_levels,
                                                      thrust::raw_pointer_cast(d_circulations_indices_no_zero.data()),
                                                      thrust::raw_pointer_cast(d_n_circulations_all_loc_ml.data()),
                                                      thrust::raw_pointer_cast(d_ce_all_location_from_indices.data()),
                                                      thrust::raw_pointer_cast(d_ce_all_location_target_indices.data()),
                                                      thrust::raw_pointer_cast(d_ce_all_moving_levels.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());
//    ThrustTVectorHost<int> h_all_location_from_indices = d_all_location_from_indices;
//    ThrustTVectorHost<int> h_all_location_target_indices = d_all_location_target_indices;
//    ThrustTVectorHost<int> h_all_moving_levels = d_all_moving_levels;
//    for(int i = 0; i < h_n_circulations_all_loc_ml.size(); i++){
//        int circulate_loc_index = i / Model::CONFIG->circulation_info().number_of_moving_levels;
//        int from_location = thrust::get<0>(h_circulations_indices_no_zero[circulate_loc_index]);
//        int target_location = thrust::get<1>(h_circulations_indices_no_zero[circulate_loc_index]);
//        int moving_level = i % Model::CONFIG->circulation_info().number_of_moving_levels;
//        int from_location2 = h_all_location_from_indices[i];
//        int target_location2 = h_all_location_target_indices[i];
//        int moving_level2 = h_all_moving_levels[i];
//        if(h_n_circulations_all_loc_ml[i] == 0) continue;
//        printf("%d from %d to %d moving level %d (%d to %d moving level %d) size %d\n",i,from_location,target_location,moving_level,
//               from_location2,target_location2,moving_level2,h_n_circulations_all_loc_ml[i]);
//    }
//    printf("\n");
    /*
     * Remove zero values in d_n_circulations_all_loc_ml
     * scan d_n_circulations_all_loc_ml and copy non-zero from & to locations from 2 vectors
     * to d_circulations_all_loc_ml_indices_no_zero
     */
    no_zero_size = thrust::count_if(d_n_circulations_all_loc_ml.begin(), d_n_circulations_all_loc_ml.end(), notZero<unsigned int>());
    auto loc_index_begin_2 = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.begin(),
                                                                          d_ce_all_location_target_indices.begin(),
                                                                          d_ce_all_moving_levels.begin()));
    auto loc_index_end_2 = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.end(),
                                                                        d_ce_all_location_target_indices.end(),
                                                                        d_ce_all_moving_levels.end()));
    ThrustT3TupleVectorDevice<int,int,int> d_circulations_all_loc_ml_indices_no_zero(no_zero_size);
    thrust::copy_if(thrust::device,
                    loc_index_begin_2,
                    loc_index_end_2,
                    d_n_circulations_all_loc_ml.begin(),
                    d_circulations_all_loc_ml_indices_no_zero.begin(),
                    copyNotZero());
    nend = thrust::remove(thrust::device,d_n_circulations_all_loc_ml.begin(),d_n_circulations_all_loc_ml.end(),0);
    ThrustTVectorDevice<unsigned int> d_n_circulations_all_loc_ml_no_zero(no_zero_size);
    thrust::copy(d_n_circulations_all_loc_ml.begin(),nend,d_n_circulations_all_loc_ml_no_zero.begin());

    /*
     * d_circulations_all_loc_ml_indices_no_zero is tuple of from_location, target_location and moving level
     * d_n_circulations_all_loc_ml_no_zero is non-zero circulation number in all locations and moving levels
     * d_n_circulations_all_loc_ml_no_zero and d_circulations_all_loc_ml_indices_no_zero are the same size and order
     * d_circulate_all_loc_ml_today is tuple of from_location, target_location, moving_level, person_index
     * */
    int total_circulations = thrust::reduce(thrust::device,d_n_circulations_all_loc_ml_no_zero.begin(),d_n_circulations_all_loc_ml_no_zero.end());
    ThrustT4TupleVectorDevice<int,int,int,unsigned int> d_circulate_all_loc_ml_today(d_circulations_all_loc_ml_indices_no_zero.size(),thrust::make_tuple(-1,-1,-1,0));
    block_size = (d_circulate_all_loc_ml_today.size() + n_threads - 1)/n_threads;
    zip_location_indices_and_n_circulations<<<block_size,n_threads>>>(d_circulate_all_loc_ml_today.size(),
                                                                      thrust::raw_pointer_cast(d_circulations_all_loc_ml_indices_no_zero.data()),
                                                                      thrust::raw_pointer_cast(d_n_circulations_all_loc_ml_no_zero.data()),
                                                                      thrust::raw_pointer_cast(d_circulate_all_loc_ml_today.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());

    thrust::sort(thrust::device,d_circulate_all_loc_ml_today.begin(),d_circulate_all_loc_ml_today.end(),circulateLess());

//    ThrustT4TupleVectorHost<int,int,int,unsigned int> h_circulate_all_loc_ml_today = d_circulate_all_loc_ml_today;
//    for(int i = 0; i < h_circulate_all_loc_ml_today.size(); i++){
//        int from_location = thrust::get<0>(h_circulate_all_loc_ml_today[i]);
//        int target_location = thrust::get<1>(h_circulate_all_loc_ml_today[i]);
//        int moving_level = thrust::get<2>(h_circulate_all_loc_ml_today[i]);
//        unsigned int size = thrust::get<3>(h_circulate_all_loc_ml_today[i]);
//        printf("%d from %d to %d moving level %d size %d\n",i,from_location,target_location,moving_level,size);
//    }
//    printf("\n");

    /*
     * Random persons based on d_n_circulations_all_loc_ml_no_zero
     * This needs to be done in 2 passes, 1st pass to fill all index with 1 person first, 1st pass is done on GPU
     * 2nd pass to fill all index with n_persons > 1, 2nd pass is done on CPU
     * In this is the first pass, fill all index in d_circulate_person_indices_today
     * with 1 person first
     * */
    ThrustT5TupleVectorDevice<int,int,int,unsigned int,int> d_circulate_person_indices_today(total_circulations,thrust::make_tuple(-1,-1,-1,-1,-1));
    block_size = (d_n_circulations_all_loc_ml_no_zero.size() + n_threads - 1)/n_threads;
    fill_circulate_person_indices<<<block_size,n_threads>>>(d_circulate_all_loc_ml_today.size(),
                                                           Model::CONFIG->circulation_info().number_of_moving_levels,
                                                           Model::GPU_RANDOM->d_states,
                                                           thrust::raw_pointer_cast(d_circulate_all_loc_ml_today.data()),
                                                           thrust::raw_pointer_cast(d_ce_popsize_by_moving_level.data()),
                                                           thrust::raw_pointer_cast(d_circulate_person_indices_today.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());

    /*
     * Check if there is any index with n_persons > 1 and random on CPU
     * otherwise schedule person events
     * */
    ThrustT5TupleVectorHost<int,int,int,unsigned int,int> h_circulate_person_indices_today = d_circulate_person_indices_today;
    auto *pi = Model::POPULATION->get_person_index<PersonIndexByLocationMovingLevel>();
    for(int i = 0; i < d_circulate_all_loc_ml_today.size();i++){
        int from_location = h_circulate_person_indices_today[i].get<0>();
        int target_location = h_circulate_person_indices_today[i].get<1>();
        int moving_level = h_circulate_person_indices_today[i].get<2>();
        int n_persons = h_circulate_person_indices_today[i].get<3>();
        auto size = static_cast<int>(pi->vPerson()[from_location][moving_level].size());
        if (size==0) continue;
        if(n_persons == 1){
            int p_index = h_circulate_person_indices_today[i].get<4>();
            Person* p = pi->vPerson()[from_location][moving_level][p_index];
            assert(p->host_state()!=Person::DEAD);
            p->today_target_locations()->push_back(target_location);
            p->randomly_choose_target_location();
//            printf("i %d GPU from %d to %d moving level %d n_persons %d p_index %d\n",
//                   i,
//                   from_location,target_location,moving_level,
//                   n_persons,p_index);
        }
        else{
            for(int j = 0; j < n_persons; j++) {
                int p_index = Model::RANDOM->random_uniform(size);
                Person* p = pi->vPerson()[from_location][moving_level][p_index];
                assert(p->host_state()!=Person::DEAD);
                p->today_target_locations()->push_back(target_location);
                p->randomly_choose_target_location();
//                printf("i %d j %d CPU from %d to %d moving level %d n_persons %d p_index %d\n",
//                       i,j,
//                       from_location,target_location,moving_level,
//                       n_persons,p_index);
            }
        }
    }

    if(Model::CONFIG->debug_config().enable_debug_text){
        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        LOG_IF(Model::SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << "[GPU Population] Update population circulation GPU (" << d_circulations_indices_no_zero.size() << " " << d_num_leavers_from_target_no_zero.size()
        << " " << total_leavers << " " << total_circulations << ") event time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms ";
    }
}

void GPU::Population::perform_infection_event() {
    auto tp_start = std::chrono::high_resolution_clock::now();
    auto tracking_index = Model::SCHEDULER->current_time() % Model::CONFIG->number_of_tracking_days();

    /*
     * Calculate probability of leaving location in all locations (n_location*n_location)
     * Also get indices from and to arrays
     * */

    int n_locations = 5;
    int n_person_each_location = 20;
    int n_sample_each_location = 4;
    TVector<double> distribution(n_locations*n_person_each_location,0.0);
    TVector<double> sum_distribution(n_locations);
    std::vector<Person*> all_persons;
    for (int i = 0; i < n_locations*n_person_each_location; ++i) {
        int location_index = i / n_person_each_location;
        auto *p = new Person();
        p->set_id(i);
        distribution[location_index*n_person_each_location+0] = 1.0;
        distribution[location_index*n_person_each_location+5] = 0.2;
        distribution[location_index*n_person_each_location+10] = 0.5;
        distribution[location_index*n_person_each_location+15] = 0.8;
        sum_distribution[location_index] = distribution[location_index*n_person_each_location+0]
                                            + distribution[location_index*n_person_each_location+5]
                                            + distribution[location_index*n_person_each_location+10]
                                            + distribution[location_index*n_person_each_location+15];
        all_persons.push_back(p);
//        printf("location %d person %d\n",location_index,p->id());
    }
    for(int i = 0; i < n_locations; i++){
        printf("location %d sum %f\n",i,sum_distribution[i]);
    }
    ThrustTVectorHost<double> h_distribution(distribution);
    ThrustTVectorHost<double> h_sum_distribution(sum_distribution);
    ThrustTVectorDevice<double> d_distribution = h_distribution;
    ThrustTVectorDevice<double> d_sum_distribution = h_sum_distribution;

    auto result = Model::RANDOM->multinomial_sampling(n_sample_each_location,
                                                   distribution,
                                                   all_persons,
                                                   false,
                                                   sum_distribution[0]);

    for(auto s : result){
        printf("CPU location %d, multinomial samples %d\n", 0, s->id());
    }

    result = Model::RANDOM->roulette_sampling(n_sample_each_location,
                                                 distribution,
                                                 all_persons,
                                                 false,
                                                 sum_distribution[1]);

    for(auto s : result){
        printf("CPU location %d, roulette samples %d\n", 1, s->id());
    }

    auto result_gpu = Model::GPU_RANDOM->roulette_sampling(n_locations,
                                                           n_sample_each_location,
                                                           all_persons,
                                                           d_distribution,
                                                           d_sum_distribution,
                                                           false);

    for(int i = 0; i <result_gpu.size(); i++){
        int location_index = i / n_sample_each_location;
        printf("GPU Roulette location %d result id %d\n",location_index,result_gpu[i]->id());
    }

    result_gpu = Model::GPU_RANDOM->multinomial_sampling(n_locations,
                                                      n_sample_each_location,
                                                      all_persons,
                                                      d_distribution,
                                                      d_sum_distribution,
                                                      false);

    for(int i = 0; i <result_gpu.size(); i++){
        int location_index = i / n_sample_each_location;
        printf("GPU Multinomial location %d result id %d\n",location_index,result_gpu[i]->id());
    }

    exit(0);
}
