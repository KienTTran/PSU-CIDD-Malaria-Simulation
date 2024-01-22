//
// Created by kient on 12/31/2023.
//

#include <curand_kernel.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/sort.h>
#include "PopulationKernel.cuh"
#include "Population.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Spatial/SpatialModel.hxx"
#include "Core/Random.h"
#include "Gpu/Utils/Utils.cuh"
#include "Gpu/Core/Random.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationMovingLevel.cuh"
#include "Properties/PersonIndexGPU.cuh"
#include "ClonalParasitePopulation.cuh"
#include "SingleHostClonalParasitePopulations.cuh"
#include "ImmuneSystem.cuh"

GPU::PopulationKernel::PopulationKernel() {
}

void GPU::PopulationKernel::init() {
    h_ie_foi_N_days_all_locations = TVector<double>(Model::CONFIG->number_of_locations()
            *Model::CONFIG->number_of_tracking_days());
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

void GPU::PopulationKernel::calculate_circulate_locations(int n_locations,ThrustTVectorDevice<int> &d_n_circulations_from_locations,ThrustTVectorDevice<double> &d_relative_outmovement_from_target_,
                                         ThrustTVectorDevice<int> &d_all_location_from_indices,ThrustTVectorDevice<int> &d_all_location_target_indices) {
    //Has to get pointer to device otherwise it will copy vector from host to device
    d_ce_popsize_residence_by_location = Model::GPU_DATA_COLLECTOR->popsize_residence_by_location();
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

void GPU::PopulationKernel::calculate_moving_level_density(ThrustTuple2VectorDevice<int,int> d_circulation_indices,ThrustTVectorDevice<double> &d_moving_level_density) {
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
void GPU::PopulationKernel::perform_circulation_event() {
    auto tp_start = std::chrono::high_resolution_clock::now();

    /*
     * Calculate probability of leaving location in all locations (n_location*n_location)
     * Also get indices from and to arrays
     * */
    Model::GPU_POPULATION_KERNEL->calculate_circulate_locations(Model::CONFIG->number_of_locations(),
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
    ThrustTuple2VectorDevice<int,int> d_circulations_indices_no_zero(no_zero_size);

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
    Model::GPU_POPULATION_KERNEL->calculate_moving_level_density(d_circulations_indices_no_zero,d_moving_level_density);

    ThrustTVectorDevice<unsigned int> d_n_circulations_all_loc_ml(d_moving_level_density.size(),0);
    Model::GPU_RANDOM->random_multinomial(d_circulations_indices_no_zero.size(),
                                         Model::CONFIG->circulation_info().number_of_moving_levels,
                                         d_num_leavers_from_target_no_zero,
                                         d_moving_level_density,
                                         d_n_circulations_all_loc_ml);
    /*
     * d_circulations_indices_no_zero and d_n_circulations_all_loc_ml not the same size
     */
//    TVector<unsigned int> h_n_circulations_all_loc_ml = d_n_circulations_all_loc_ml;
//    ThrustT2TupleVector<int,int> h_circulations_indices_no_zero = d_circulations_indices_no_zero;
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
//    TVector<int> h_all_location_from_indices = d_all_location_from_indices;
//    TVector<int> h_all_location_target_indices = d_all_location_target_indices;
//    TVector<int> h_all_moving_levels = d_all_moving_levels;
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
    ThrustTuple3VectorDevice<int,int,int> d_circulations_all_loc_ml_indices_no_zero(no_zero_size);
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
    ThrustTuple4VectorDevice<int,int,int,unsigned int> d_circulate_all_loc_ml_today(d_circulations_all_loc_ml_indices_no_zero.size(),thrust::make_tuple(-1,-1,-1,0));
    block_size = (d_circulate_all_loc_ml_today.size() + n_threads - 1)/n_threads;
    zip_location_indices_and_n_circulations<<<block_size,n_threads>>>(d_circulate_all_loc_ml_today.size(),
                                                                      thrust::raw_pointer_cast(d_circulations_all_loc_ml_indices_no_zero.data()),
                                                                      thrust::raw_pointer_cast(d_n_circulations_all_loc_ml_no_zero.data()),
                                                                      thrust::raw_pointer_cast(d_circulate_all_loc_ml_today.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());

    thrust::sort(thrust::device,d_circulate_all_loc_ml_today.begin(),d_circulate_all_loc_ml_today.end(),circulateLess());

//    ThrustT4TupleVector<int,int,int,unsigned int> h_circulate_all_loc_ml_today = d_circulate_all_loc_ml_today;
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
    ThrustTuple5VectorDevice<int,int,int,unsigned int,int> d_circulate_person_indices_today(total_circulations,thrust::make_tuple(-1,-1,-1,-1,-1));
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
    ThrustTuple5VectorHost<int,int,int,unsigned int,int> h_circulate_person_indices_today = d_circulate_person_indices_today;
    auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationMovingLevel>();
    for(int i = 0; i < d_circulate_all_loc_ml_today.size();i++){
        int from_location = h_circulate_person_indices_today[i].get<0>();
        int target_location = h_circulate_person_indices_today[i].get<1>();
        int moving_level = h_circulate_person_indices_today[i].get<2>();
        int n_persons = h_circulate_person_indices_today[i].get<3>();
        auto size = static_cast<int>(pi->vPerson()[from_location][moving_level].size());
        if (size==0) continue;
        if(n_persons == 1){
            int p_index = h_circulate_person_indices_today[i].get<4>();
            GPU::Person* p = pi->vPerson()[from_location][moving_level][p_index];
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
                GPU::Person* p = pi->vPerson()[from_location][moving_level][p_index];
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
        LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << "[GPU Population] Update population circulation GPU (" << d_circulations_indices_no_zero.size() << " " << d_num_leavers_from_target_no_zero.size()
        << " " << total_leavers << " " << total_circulations << ") event time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms ";
    }
}

void GPU::PopulationKernel::calculate_n_person_bitten_today(int n_locations,
                                                      ThrustTVectorDevice<double> &d_foi_all_locations,
                                                      ThrustTVectorDevice<int> &d_n_person_bitten_today_all_locations){


}

void GPU::PopulationKernel::perform_infection_event() {
    auto tp_start = std::chrono::high_resolution_clock::now();
    auto tracking_index = Model::GPU_SCHEDULER->current_time() % Model::CONFIG->number_of_tracking_days();

    /*
     * Calculate probability of leaving location in all locations (n_location*n_location)
     * Also get indices from and to arrays
     * */

    ThrustTVectorDevice<double> d_foi_all_locations;
    ThrustTVectorDevice<int> d_n_person_bitten_today_all_locations;
    calculate_n_person_bitten_today(Model::CONFIG->number_of_locations(),
                                    d_foi_all_locations,d_n_person_bitten_today_all_locations);
}

__global__ void update_current_foi_kernel(int size){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < size; index += stride) {
    }
}

void GPU::PopulationKernel::update_current_foi(){

}

/*
 * In order to use virtual in kernel, the base class needs to be instanced on GPU
 * */
__global__ void set_gpu_update_function_kernel(int size,
                                               GPU::ParasiteDensityUpdateFunction** d_function,
//                                               GPU::ImmuneSystem** d_immune_system,
                                               int* type){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < size; index += stride) {
        switch(type[index]) {
            case 1:
                d_function[index] = new GPU::ClinicalUpdateFunction();
                break;
            case 2:
                d_function[index] = new GPU::ImmunityClearanceUpdateFunction();
                break;
            default:
                printf("[GPU::ClonalParasitePopulation] ERROR: GPU::ParasiteDensityUpdateFunction not set\n");
                break;
        }
//        d_immune_system[index] = new GPU::ImmuneSystem();
    }
}

/*
 * Set update function inside kernel
 * To do this the base class of update function needs to be on GPU first
 * Remember to run
 * h_cpp->set_h_parasite_density_level(Model::CONFIG->parasite_density_level());
 * h_cpp->allocate_on_gpu();
 * before calling this function
 * */
//__global__ void update_all_individuals_kernel(int size,
//                                              int current_time,
//                                              ParasiteDensityLevel h_parasite_density_level,
//                                              ImmuneSystemInformation h_immune_system_information,
////                                              GPU::ParasiteDensityUpdateFunction** d_update_function,
////                                              GPU::ImmuneSystem** d_immune_system,
//                                              GPU::ClonalParasitePopulation** d_cpp){
//    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
//    int stride = blockDim.x * gridDim.x;
//    for(int index = thread_index; index < size; index += stride) {
////        printf("index %d cpp %d\n",index,d_update_function[index]->type());
//        d_cpp[index]->update_gpu(h_parasite_density_level,h_immune_system_information,
////                          d_update_function[index],d_immune_system[index],
//                          current_time);
//    }
//}


/*
 * Set update function inside kernel
 * To do this the base class of update function needs to be on GPU first
 * Remember to run
 * h_cpp->set_h_parasite_density_level(Model::CONFIG->parasite_density_level());
 * h_cpp->allocate_on_gpu();
 * before calling this function
 * */
__global__ void update_all_individuals_kernel2(int size,
                                              int* d_person_id,
                                              int current_time,
                                              int* lastest_updated_time,
                                              ParasiteDensityLevel h_parasite_density_level,
                                              ImmuneSystemInformation *h_immune_system_information,
                                              GPU::ClonalParasitePopulation** d_cpp,
                                              GPU::Genotype** d_genotype,
                                              GPU::ParasiteDensityUpdateFunction** d_update_function,
                                              GPU::ImmuneSystem** d_immune_system){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < size; index += stride) {
        if (lastest_updated_time[index] == current_time) return;
        d_person_id[index],d_cpp[index]->get_current_parasite_density_gpu(d_update_function[index],
                                                                          d_genotype[index],
                                                                          d_immune_system[index],
                                                                          h_parasite_density_level,
                                                                          h_immune_system_information,
                                                                          current_time,
                                                                          lastest_updated_time[index]);
//        printf("kernel d_cpp %d %f\n",d_person_id[index],d_cpp[index]->get_current_parasite_density_gpu(d_update_function[index],
//                                                                                                        d_genotype[index],
//                                                                                                        d_immune_system[index],
//                                                                                                        h_parasite_density_level,
//                                                                                                        h_immune_system_information,
//                                                                                                        current_time,
//                                                                                                        lastest_updated_time[index]));
        lastest_updated_time[index] = current_time;
    }
}

void GPU::PopulationKernel::update_all_individuals(){
    auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
    TVector<GPU::ClonalParasitePopulation*> h_cpp;
    TVector<GPU::Genotype*> h_genotype;
    TVector<GPU::ImmuneSystem*> h_immune_system;
    TVector<int> h_cpp_update_function_type;
    TVector<int> h_person_id;
    TVector<int> h_latest_update_time;
    for(int i = 0; i < pi->h_persons().size(); i++){
        auto sh = pi->h_persons()[i]->all_clonal_parasite_populations();
        if(sh->size() == 0) continue;
        for(int j = 0; j < sh->size(); j++){
            if(sh->parasites()->at(j)->update_function() == nullptr) continue;
            if(sh->parasites()->at(j)->update_function()->type() == 0) continue;
            h_genotype.push_back(sh->parasites()->at(j)->genotype());
            h_immune_system.push_back(pi->h_persons()[i]->immune_system());
            h_cpp.push_back(sh->parasites()->at(j));
            h_cpp_update_function_type.push_back(sh->parasites()->at(j)->update_function()->type());
            h_person_id.push_back(pi->h_persons()[i]->id());
            h_latest_update_time.push_back(pi->h_persons()[i]->latest_update_time());
        }
    }
    if(h_cpp.size() == 0) return;

    printf("h_cpp size %d\n",h_cpp.size());
    const int h_cpp_size = h_cpp.size();

//    printf("Host:\n");
//    for (int i = 0; i < h_cpp_size; i++) {
//        std::cout << h_genotype[i]->test2() << std::endl;
//        std::cout << h_immune_system[i]->test_ << std::endl;
//        printf("%d type %d\n",i,h_cpp_update_function_type[i]);
//    }

    /*
     * Since ParasiteDensityUpdateFunction and ImmuneSystem are virtual classes,
     * they need to be instanced on GPU first
     * */
    int *d_cpp_update_function_type;
    cudaMalloc((void**)&d_cpp_update_function_type, sizeof(int) * h_cpp_size);
    cudaMemcpy(d_cpp_update_function_type, h_cpp_update_function_type.data(), sizeof(int) * h_cpp_size, cudaMemcpyHostToDevice);
    check_cuda_error(cudaGetLastError());

    GPU::ParasiteDensityUpdateFunction **d_update_function;
    cudaMalloc((void**)&d_update_function, h_cpp_size*sizeof(GPU::ParasiteDensityUpdateFunction*));
    check_cuda_error(cudaGetLastError());

//    GPU::ImmuneSystem **d_immune_system;
//    cudaMalloc((void**)&d_immune_system, h_cpp_size*sizeof(GPU::ImmuneSystem*));
//    check_cuda_error(cudaGetLastError());

    int n_threads = 256;
    int block_size = (h_cpp_size + n_threads - 1)/n_threads;
    set_gpu_update_function_kernel<<<block_size,n_threads>>>(h_cpp_size,
                                            d_update_function,
//                                            d_immune_system,
                                            d_cpp_update_function_type);
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());


//    /*
//     * ClonalParasitePopulation and Genotype are not virtual classes,
//     * they can be copied from CPU
//     * refer to https://stackoverflow.com/questions/38978297/pointer-to-array-of-pointers-to-objects-in-cuda
//     * */
//    GPU::ClonalParasitePopulation** h_d_cpp = (GPU::ClonalParasitePopulation**)malloc(sizeof(GPU::ClonalParasitePopulation*) * h_cpp_size);
//    for (int i = 0; i < h_cpp_size; i++) {
//        // Allocate space for an Obj and assign
//        cudaMalloc((void**)&h_d_cpp[i], sizeof(GPU::ClonalParasitePopulation));
//        // Copy the object to the device (only has single scalar field to keep it simple)
//        cudaMemcpy(h_d_cpp[i], &(h_cpp[i]), sizeof(GPU::ClonalParasitePopulation), cudaMemcpyHostToDevice);
//        check_cuda_error(cudaGetLastError());
//    }
//
//    // Create a pointer which will point to device memory
//    GPU::ClonalParasitePopulation** d_d_cpp = NULL;
//    // Allocate space for 3 pointers on device at above location
//    cudaMalloc((void**)&d_d_cpp, sizeof(GPU::ClonalParasitePopulation*) * h_cpp_size);
//    // Copy the pointers from the host memory to the device array
//    cudaMemcpy(d_d_cpp, h_d_cpp, sizeof(GPU::ClonalParasitePopulation*) * h_cpp_size, cudaMemcpyHostToDevice);
//    check_cuda_error(cudaGetLastError());
//
//    GPU::Genotype** h_d_genotype = (GPU::Genotype**)malloc(sizeof(GPU::Genotype*) * h_cpp_size);
//    for (int i = 0; i < h_cpp_size; i++) {
//        // Allocate space for an Obj and assign
//        cudaMalloc((void**)&h_d_genotype[i], sizeof(GPU::Genotype));
//        // Copy the object to the device (only has single scalar field to keep it simple)
//        cudaMemcpy(h_d_genotype[i], &(h_genotype[i]), sizeof(GPU::Genotype), cudaMemcpyHostToDevice);
//        check_cuda_error(cudaGetLastError());
//    }
//
//    // Create a pointer which will point to device memory
//    GPU::Genotype** d_d_genotype = NULL;
//    // Allocate space for 3 pointers on device at above location
//    cudaMalloc((void**)&d_d_genotype, sizeof(GPU::Genotype*) * h_cpp_size);
//    // Copy the pointers from the host memory to the device array
//    cudaMemcpy(d_d_genotype, h_d_genotype, sizeof(GPU::Genotype*) * h_cpp_size, cudaMemcpyHostToDevice);
//    check_cuda_error(cudaGetLastError());
//
//    ImmuneSystemInformation *d_immune_system_information;
//    cudaMalloc((void**)&d_immune_system_information, sizeof(ImmuneSystemInformation));
//    cudaMemcpy(d_immune_system_information, &Model::CONFIG->immune_system_information(), sizeof(ImmuneSystemInformation), cudaMemcpyHostToDevice);
//    check_cuda_error(cudaGetLastError());
//
//    int *d_person_id;
//    cudaMalloc((void**)&d_person_id, sizeof(int) * h_cpp_size);
//    cudaMemcpy(d_person_id, h_person_id.data(), sizeof(int) * h_cpp_size, cudaMemcpyHostToDevice);
//    check_cuda_error(cudaGetLastError());
//
//    int *d_latest_update_time;
//    cudaMalloc((void**)&d_latest_update_time, sizeof(int) * h_cpp_size);
//    cudaMemcpy(d_latest_update_time, h_latest_update_time.data(), sizeof(int) * h_cpp_size, cudaMemcpyHostToDevice);
//    check_cuda_error(cudaGetLastError());
//
////    update_all_individuals_kernel2<<<block_size,n_threads>>>(h_cpp_size,
////                                                            d_person_id,
////                                                            Model::GPU_SCHEDULER->current_time(),
////                                                            d_latest_update_time,
////                                                            Model::CONFIG->parasite_density_level(),
////                                                            d_immune_system_information,
////                                                            d_d_cpp,
////                                                            d_d_genotype,
////                                                            d_update_function,
////                                                            d_immune_system);
////    cudaDeviceSynchronize();
////    check_cuda_error(cudaGetLastError());
//
////    for (int i = 0; i < h_cpp_size; i++) {
////        cudaMemcpy(h_cpp[i], h_d_cpp[i], sizeof(GPU::ClonalParasitePopulation), cudaMemcpyDeviceToHost);
//////        cudaMemcpy(h_genotype[i], h_d_genotype[i], sizeof(GPU::Genotype), cudaMemcpyDeviceToHost);
////        cudaMemcpy(h_immune_system[i], h_d_immune_system[i], sizeof(GPU::ImmuneSystem), cudaMemcpyDeviceToHost);
////        check_cuda_error(cudaGetLastError());
////    }
//
//    // Write out
////    printf("D2H Host:\n");
////    for (int i = 0; i < h_cpp_size; i++) {
////        std::cout << h_cpp[i]->test2() << std::endl;
////        std::cout << h_genotype[i]->test_ << std::endl;
////    }
//
//    for(int i = 0; i < h_cpp_size; i++){
//        cudaFree(h_d_cpp[i]);
//        cudaFree(h_d_genotype[i]);
//    }
//    cudaFree(d_immune_system);
//    cudaFree(d_update_function);
//    cudaFree(d_immune_system_information);
//    cudaFree(d_person_id);
//    cudaFree(d_latest_update_time);
//    h_cpp.clear();
//    h_genotype.clear();
//    h_immune_system.clear();
}
