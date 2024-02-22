//
// Created by kient on 12/31/2023.
//

#include <curand_kernel.h>
#include <cuda_profiler_api.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/sort.h>
#include <thrust/iterator/discard_iterator.h>
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
#include "Gpu/Therapies/DrugDatabase.cuh"
#include <math.h>

#ifndef DRUG_CUT_OFF_VALUE
#define DRUG_CUT_OFF_VALUE 0.1
#endif

/*
 * From Random.cu and Utils.cu
 * */

__device__ double curand_gsl_uniform_double(double rand, double min, double max);

__device__ int curand_gsl_uniform_int(double rand, int min, int max);

__device__ double curand_gsl_cdf_ugaussian_P(double x);

__device__ bool is_equal(double a, double b, double epsilon);

__host__ int encode_vec2_to_int(TVector<int> int_vector);

__device__ int *decode_int_to_arr2(int encoded_value, int *result);

GPU::PopulationKernel::PopulationKernel() {
}

void GPU::PopulationKernel::init() {
  /*
   * To access resistant_aa_locations in device, we need to copy the resistant_aa_locations from host to device
   * and make an index map to access the resistant_aa_locations. h_drug_res_aa_loc_index contain start index and
   * end index of res_aa_loc of each drug.
   * */
  h_drug_res_aa_loc_index = TVector<int>(Model::CONFIG->gpu_drug_db()->size(), -1);
  int start_index = 0;
  int count = 0;
  for (int d_index = 0; d_index < Model::CONFIG->gpu_drug_db()->size(); d_index++) {
    if (Model::CONFIG->gpu_drug_db()->at(d_index)->resistant_aa_locations.size() > 0) {
      for (auto res_aa: Model::CONFIG->gpu_drug_db()->at(d_index)->resistant_aa_locations) {
        h_drug_res_aa_loc.push_back(res_aa);
        printf("drug_index %d start_index: %d count: %d res_aa_chr_index: %d res_aa_gene_index: %d res_aa_aa_index: %d res_aa_index_in_string: %d\n",
               d_index,
               start_index,
               count,
               res_aa.chromosome_id,
               res_aa.gene_id,
               res_aa.aa_id,
               res_aa.aa_index_in_aa_string
        );
        count++;
      }
      h_drug_res_aa_loc_index[d_index] = encode_vec2_to_int(TVector<int>{start_index,
                                                                         static_cast<int>(start_index + Model::CONFIG->gpu_drug_db()->at(
                                                                             d_index)->resistant_aa_locations.size())});
      printf("drug_index %d (%d -> %d) -> %d\n",
             d_index,
             start_index,
             start_index + Model::CONFIG->gpu_drug_db()->at(d_index)->resistant_aa_locations.size(),
             h_drug_res_aa_loc_index[d_index]);
      start_index += Model::CONFIG->gpu_drug_db()->at(d_index)->resistant_aa_locations.size();
    }
  }
  d_drug_res_aa_loc = h_drug_res_aa_loc;
  d_drug_res_aa_loc_index = h_drug_res_aa_loc_index;
  check_cuda_error(cudaMalloc(&d_gen_mutation_mask, sizeof(char) * Model::CONFIG->mutation_mask().size()));
  check_cuda_error(cudaMemcpy(d_gen_mutation_mask, Model::CONFIG->mutation_mask().c_str(),
                              sizeof(char) * Model::CONFIG->mutation_mask().size(), cudaMemcpyHostToDevice));

  /*
   * To access genotype_max_copies in device, we need to copy the genotype_max_copies from host to device
   * and make an index map to access the genotype_max_copies
   * To access genotype_aa in device, we need to copy the genotype_aa from host to device. h_gen_aa_int_start_index
   * contains start index and end index of genotype_aa of each chromosome.
   *
   * */
  start_index = 0;
  std::string gen_aa_test;
  TVector<int> gen_aa;
  TVector<int> gen_aa_size;
  TVector<int> gen_max_copies;
  h_gen_aa_size = TVector<int>(Model::CONFIG->pf_genotype_info().chromosome_infos.size(), -1);
  h_gen_gene_size = TVector<int>(Model::CONFIG->pf_genotype_info().chromosome_infos.size(), -1);
  h_gen_max_copies = TVector<int>(Model::CONFIG->pf_genotype_info().chromosome_infos.size(), -1);
  h_gen_aa_int_start_index = TVector<int>(Model::CONFIG->pf_genotype_info().chromosome_infos.size(), -1);
  for (int chr_index = 0; chr_index < Model::CONFIG->pf_genotype_info().chromosome_infos.size(); chr_index++) {
    if (Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index].gene_infos.size() > 0) {
      h_gen_gene_size[chr_index] = Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index].gene_infos.size();
      gen_aa_size.clear();
      gen_max_copies.clear();
      int gen_aa_per_gene_size = 0;
      for (int gen_index = 0; gen_index < Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index].gene_infos.size(); gen_index++) {
        for (int aa_pos_index = 0; aa_pos_index < Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
            .gene_infos[gen_index].aa_position_infos.size(); aa_pos_index++) {
          gen_aa.clear();
          gen_aa_test = "";
          for (int aa_index = 0; aa_index < Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
              .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids.size(); aa_index++) {
            /*
             * Char to Int
             * */
            gen_aa.push_back((int) (Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
                                        .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids[aa_index] - 48));
          }
          std::string a = "";
          a = Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
              .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids[0];
          std::string b = "";
          b = Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
              .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids[1];
          std::string c = "";
          c = std::to_string((int) (Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
                                        .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids[0] - 48));
          std::string d = "";
          d = std::to_string((int) (Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index]
                                        .gene_infos[gen_index].aa_position_infos[aa_pos_index].amino_acids[1] - 48));
          gen_aa_test = a + "->" + b + " " + c + "->" + d;
          if (gen_aa.size() > 0) {
            printf("gen_aa_int: %s %d\n", gen_aa_test.c_str(), encode_vec2_to_int(gen_aa));
            h_gen_aa_int.push_back(encode_vec2_to_int(gen_aa));
            gen_aa_per_gene_size++;
          }
        }
        gen_aa_size.push_back(Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index].gene_infos[gen_index].aa_position_infos.size());
        gen_max_copies.push_back(Model::CONFIG->pf_genotype_info().chromosome_infos[chr_index].gene_infos[gen_index].max_copies);
      }
      /*
       * This is just for 2 gene in 1 chromosome, if number of gene is bigger then 2 then we need to find another way
       * */
      if (gen_aa_size.size() > 1) {
        h_gen_aa_size[chr_index] = encode_vec2_to_int(gen_aa_size);
        h_gen_max_copies[chr_index] = encode_vec2_to_int(gen_max_copies);
      } else {
        h_gen_aa_size[chr_index] = gen_aa_size[0];
        h_gen_max_copies[chr_index] = gen_max_copies[0];
      }
      if (h_gen_aa_size[chr_index] > 0) {
        h_gen_aa_int_start_index[chr_index] = start_index;
        start_index += gen_aa_per_gene_size;
      }
    }
  }
  d_gen_gene_size = h_gen_gene_size;
  d_gen_aa_size = h_gen_aa_size;
  d_gen_max_copies = h_gen_max_copies;
  d_gen_aa_int = h_gen_aa_int;
  d_gen_aa_int_start_index = h_gen_aa_int_start_index;

  /*
   * Here we don't copy Model::CONFIG->parasite_density_level()
   * to device because there is no vector or pointers in this struct
   * If it has any vector or pointers, we need to copy it to device first.
   * In the same manner we will need to copy ImmuneSystemInformation to device.
   * */
  check_cuda_error(cudaMalloc((void **) &d_immune_system_information, sizeof(ImmuneSystemInformation)));
  check_cuda_error(cudaMemcpy(d_immune_system_information, &Model::CONFIG->immune_system_information(),
                              sizeof(ImmuneSystemInformation), cudaMemcpyHostToDevice));
  check_cuda_error(cudaGetLastError());

  pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
  d_streams = new cudaStream_t[Model::CONFIG->gpu_config().n_streams];
  for (int i = 0; i < Model::CONFIG->gpu_config().n_streams; ++i) {
    check_cuda_error(cudaStreamCreate(&d_streams[i]));
  }
  check_cuda_error(cudaEventCreate(&start_event));
  check_cuda_error(cudaEventCreate(&stop_event));

  h_sum_biting_moving_foi_by_loc = TVector<ThrustTuple4<int,double,double,double>>(Model::CONFIG->number_of_locations(),
      ThrustTuple4<int,double,double,double>(0,0.0,0.0,0.0));
  force_of_infection_for_N_days_by_location = TVector<TVector<double>>(
      Model::CONFIG->number_of_tracking_days(), TVector<double>(Model::CONFIG->number_of_locations(), 0));
  individual_relative_biting_by_location = TVector<TVector<double>>(Model::CONFIG->number_of_locations(), TVector<double>());
  individual_relative_moving_by_location = TVector<TVector<double>>(Model::CONFIG->number_of_locations(), TVector<double>());
  individual_foi_by_location = TVector<TVector<double>>(Model::CONFIG->number_of_locations(), TVector<double>());
}

/*
 * Calculate number of circulations from each location
 * this if from poisson distribution with mean = popsize_residence_by_location * circulation_percent
 * d_n_circulations_from_locations is number of circulations from each location
 * */
__global__ void calculate_circulation_number(int n_locations,
                                             int *d_popsize_residence_by_location,
                                             double circulation_percent,
                                             curandState *d_state,
                                             int *d_n_circulations_from_locations) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  curandState local_state = d_state[thread_index];
  for (int index = thread_index; index < n_locations; index += stride) {
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
                                                    int *d_popsize_residence_by_location,
                                                    double *d_spatial_model_parameters,
                                                    double *d_spatial_model_travels,
                                                    int *d_district_vector,
                                                    double *d_distance_vector,
                                                    double *d_relative_outmovement,
                                                    int *d_from_indices,
                                                    int *d_target_indices) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int index = thread_index; index < n_locations * n_locations; index += stride) {
    int from_location = index / n_locations;
    int target_location = index % n_locations;
    if (d_distance_vector[from_location * n_locations + target_location] == 0.0) return;
    double distance = d_distance_vector[from_location * n_locations + target_location];
    double kernel = pow(1 + (distance / d_spatial_model_parameters[2]), (-d_spatial_model_parameters[1]));
    double probability = pow(d_popsize_residence_by_location[from_location], d_spatial_model_parameters[0]) * kernel;
    probability = probability / (1 + d_spatial_model_travels[from_location] + d_spatial_model_travels[target_location]);
    if (d_district_vector[from_location] == static_cast<int>(d_spatial_model_parameters[3]) &&
        d_district_vector[target_location] == static_cast<int>(d_spatial_model_parameters[3])) {
      probability /= d_spatial_model_parameters[4];
    }
    d_relative_outmovement[index] = probability;
    d_from_indices[index] = from_location;
    d_target_indices[index] = target_location;
  }
}

void GPU::PopulationKernel::calculate_circulate_locations(int n_locations, ThrustTVectorDevice<int> &d_n_circulations_from_locations,
                                                          ThrustTVectorDevice<double> &d_relative_outmovement_from_target_,
                                                          ThrustTVectorDevice<int> &d_all_location_from_indices,
                                                          ThrustTVectorDevice<int> &d_all_location_target_indices) {
  //Has to get pointer to device otherwise it will copy vector from host to device
  d_ce_popsize_residence_by_location = Model::GPU_DATA_COLLECTOR->popsize_residence_by_location();
  d_ce_spatial_model_parameters = Model::CONFIG->spatial_model()->getSpatialModelParameters();
  d_ce_spatial_model_travels = Model::CONFIG->spatial_model()->getSpatialModelTravels();
  d_ce_spatial_districts = Model::CONFIG->h_spatial_districts;
  d_ce_spatial_distances = Model::CONFIG->h_spatial_distances;

  //All probabilities because thrust run all arrays at the same time. If use 1 array then values are overwritten.
  d_n_circulations_from_locations.resize(n_locations, 0);

  //Get circulations by location
  int n_threads = Model::CONFIG->gpu_config().n_threads;
  int block_size = ceil((d_n_circulations_from_locations.size() + n_threads - 1) / n_threads);
  calculate_circulation_number<<<block_size, n_threads>>>(n_locations,
                                                          thrust::raw_pointer_cast(d_ce_popsize_residence_by_location.data()),
                                                          Model::CONFIG->circulation_info().circulation_percent,
                                                          Model::GPU_RANDOM->d_states,
                                                          thrust::raw_pointer_cast(d_n_circulations_from_locations.data()));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());

  //Get outmovement probabilities
  d_relative_outmovement_from_target_.resize(n_locations * n_locations, 0.0);
  d_all_location_from_indices.resize(n_locations * n_locations, 0);
  d_all_location_target_indices.resize(n_locations * n_locations, 0);
  n_threads = Model::CONFIG->gpu_config().n_threads;
  block_size = (d_relative_outmovement_from_target_.size() + n_threads - 1) / n_threads;
  calculate_circulation_probabilities<<<block_size, n_threads>>>(Model::CONFIG->number_of_locations(),
                                                                 thrust::raw_pointer_cast(d_ce_popsize_residence_by_location.data()),
                                                                 thrust::raw_pointer_cast(d_ce_spatial_model_parameters.data()),
                                                                 thrust::raw_pointer_cast(d_ce_spatial_model_travels.data()),
                                                                 thrust::raw_pointer_cast(d_ce_spatial_districts.data()),
                                                                 thrust::raw_pointer_cast(d_ce_spatial_distances.data()),
                                                                 thrust::raw_pointer_cast(d_relative_outmovement_from_target_.data()),
                                                                 thrust::raw_pointer_cast(d_all_location_from_indices.data()),
                                                                 thrust::raw_pointer_cast(d_all_location_target_indices.data()));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());
}

/*
 * Calculate moving level density at each location, output size is n_location*n_moving_level
 * This value is for each destination locations
 * */
__global__ void calculate_moving_level_density_kernel(int n_locations,
                                                      int n_moving_levels,
                                                      thrust::tuple<int, int> *d_circulation_indices,
                                                      int *d_popsize_by_moving_level,
                                                      double *d_moving_level_value,
                                                      double *d_moving_level_density) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int index = thread_index; index < n_locations * n_moving_levels; index += stride) {
    int circulate_loc_index = index / n_moving_levels;
    int from_location = thrust::get<0>(d_circulation_indices[circulate_loc_index]);
    int moving_level = index % n_moving_levels;
    if (d_popsize_by_moving_level[from_location * n_moving_levels + moving_level] == 0 || d_moving_level_value[moving_level] == 0.0) return;
    d_moving_level_density[index] = d_popsize_by_moving_level[from_location * n_moving_levels + moving_level] * d_moving_level_value[moving_level];
  }
}

void GPU::PopulationKernel::calculate_moving_level_density(ThrustTuple2VectorDevice<int, int> d_circulation_indices,
                                                           ThrustTVectorDevice<double> &d_moving_level_density) {
  if (d_circulation_indices.size() == 0) {
    return;
  }
  d_moving_level_density.resize(d_circulation_indices.size() * Model::CONFIG->circulation_info().number_of_moving_levels);
  d_ce_popsize_by_moving_level = Model::CONFIG->h_popsize_by_moving_level;
  d_ce_moving_level_value = Model::CONFIG->circulation_info().v_moving_level_value;
  int n_threads = Model::CONFIG->gpu_config().n_threads;
  int block_size = ceil((d_moving_level_density.size() + n_threads - 1) / n_threads);
  calculate_moving_level_density_kernel<<<block_size, n_threads>>>(d_circulation_indices.size(),
                                                                   Model::CONFIG->circulation_info().number_of_moving_levels,
                                                                   thrust::raw_pointer_cast(d_circulation_indices.data()),
                                                                   thrust::raw_pointer_cast(d_ce_popsize_by_moving_level.data()),
                                                                   thrust::raw_pointer_cast(d_ce_moving_level_value.data()),
                                                                   thrust::raw_pointer_cast(d_moving_level_density.data()));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());
}

template<typename T>
struct NotZero : public thrust::unary_function<T, bool> {
    __host__ __device__
    bool operator()(T x) {
      return x != 0;
    }
};

struct CopyNotZero : public thrust::unary_function<unsigned int, bool> {
    __host__ __device__
    bool operator()(unsigned int x) {
      return x != 0;
    }
};

struct CirculateLessThan {
    __host__ __device__
    bool operator()(const thrust::tuple<int, int, int, unsigned int> &t1, thrust::tuple<int, int, int, unsigned int> &t2) {
      if (t1.get<3>() < t2.get<3>())
        return true;
      if (t1.get<3>() > t2.get<3>())
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
                                                    thrust::tuple<int, int> *d_circulation_indices,
                                                    unsigned int *d_n_circulations_all_loc_ml,
                                                    int *d_from_indices,
                                                    int *d_target_indices,
                                                    int *d_moving_levels) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int index = thread_index; index < n_locations * n_moving_levels; index += stride) {
    if (d_n_circulations_all_loc_ml[index] == 0) return;
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
                                                        thrust::tuple<int, int, int> *d_circulations_all_loc_ml_indices_no_zero,
                                                        unsigned int *d_n_circulations_all_loc_ml_no_zero,
                                                        thrust::tuple<int, int, int, unsigned int> *d_circulate_all_loc_ml_today) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int index = thread_index; index < size; index += stride) {
    if (d_n_circulations_all_loc_ml_no_zero[index] == 0) return;
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
                                              thrust::tuple<int, int, int, unsigned int> *d_circulate_all_loc_ml_today,
                                              int *d_popsize_by_loc_ml,
                                              thrust::tuple<int, int, int, unsigned int, int> *d_circulate_person_indices_today) {
  int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  curandState local_state = d_state[thread_index];
  for (int index = thread_index; index < n_locations; index += stride) {
    if (thrust::get<3>(d_circulate_all_loc_ml_today[index]) == 0) return;
    int from_location = thrust::get<0>(d_circulate_all_loc_ml_today[index]);
    int target_location = thrust::get<1>(d_circulate_all_loc_ml_today[index]);
    int moving_level = thrust::get<2>(d_circulate_all_loc_ml_today[index]);
    unsigned int n_persons = thrust::get<3>(d_circulate_all_loc_ml_today[index]);
    int size = d_popsize_by_loc_ml[from_location * n_moving_levels + moving_level];
    if (size == 0 || n_persons == 0) return;
    /*
     * Random uniform to get person index
     * To get result same as gsl, using casting method, which is [0,n-1]
     * ceiling methos is [1,n]
     * https://github.com/nglee/so_answers/blob/master/cuda/170426/kernel.cu
     * */
//        printf("kernel %d from_location: %d, moving_level: %d, popsize_loc_ml: %d size: %d\n",index,
//               from_location,moving_level,popsize_loc_ml,d_n_circulations_all_loc_ml_no_zero[index]);
    if (n_persons == 1) {
      int p_index = curand_uniform(&local_state) * size;
      d_circulate_person_indices_today[index] = thrust::make_tuple(from_location, target_location, moving_level, n_persons, p_index);
//        printf("kernel %d from_location: %d, moving_level: %d, popsize_loc_ml: %d, size: %d, p_index: %d\n",index,from_location,moving_level,
//               popsize_loc_ml,d_n_circulations_all_loc_ml_no_zero[index],p_index);
    } else {
      d_circulate_person_indices_today[index] = thrust::make_tuple(from_location, target_location, moving_level, n_persons, -1);
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
  size_t no_zero_size = thrust::count_if(d_num_leavers_from_target_.begin(), d_num_leavers_from_target_.end(), NotZero<unsigned int>());
  /*
   * Remove zero values in d_num_leavers_from_target_ and d_n_circulations_from_locations
   * d_circulations_indices_no_zero index is not location index, its index is index of n_locations*n_locations
   * thrust::get<0>(d_circulations_indices_no_zero[i]) is from location index
   * thrust::get<1>(d_circulations_indices_no_zero[i]) is to location index
   * */
  ThrustTuple2VectorDevice<int, int> d_circulations_indices_no_zero(no_zero_size);

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
                  CopyNotZero());
  thrust::device_vector<unsigned int>::iterator nend = thrust::remove(thrust::device, d_num_leavers_from_target_.begin(), d_num_leavers_from_target_.end(), 0);
  ThrustTVectorDevice<unsigned int> d_num_leavers_from_target_no_zero(no_zero_size);
  thrust::copy(d_num_leavers_from_target_.begin(), nend, d_num_leavers_from_target_no_zero.begin());

  /*
   * d_num_leavers_from_target_no_zero is non-zero leavers in all locations
   * */
  int total_leavers = thrust::reduce(thrust::device, d_num_leavers_from_target_no_zero.begin(), d_num_leavers_from_target_no_zero.end());
//    LOG_IF(total_leavers == 0, DEBUG) << "[Population] Update population circulation GPU total_leavers = 0";
  if (total_leavers == 0) return;

  /*
   * Calculate moving level density at each location, output size is n_location*n_moving_level
   * n_location is no-zero indices
   */
  ThrustTVectorDevice<double> d_moving_level_density;
  Model::GPU_POPULATION_KERNEL->calculate_moving_level_density(d_circulations_indices_no_zero, d_moving_level_density);

  ThrustTVectorDevice<unsigned int> d_n_circulations_all_loc_ml(d_moving_level_density.size(), 0);
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
  if (d_n_circulations_all_loc_ml.size() == 0) return;
  int n_threads = Model::CONFIG->gpu_config().n_threads;
  int block_size = ceil((d_n_circulations_all_loc_ml.size() + n_threads - 1) / n_threads);
  extract_locations_and_moving_levels<<<block_size, n_threads>>>(d_circulations_indices_no_zero.size(),
                                                                 Model::CONFIG->circulation_info().number_of_moving_levels,
                                                                 thrust::raw_pointer_cast(d_circulations_indices_no_zero.data()),
                                                                 thrust::raw_pointer_cast(d_n_circulations_all_loc_ml.data()),
                                                                 thrust::raw_pointer_cast(d_ce_all_location_from_indices.data()),
                                                                 thrust::raw_pointer_cast(d_ce_all_location_target_indices.data()),
                                                                 thrust::raw_pointer_cast(d_ce_all_moving_levels.data()));
  check_cuda_error(cudaDeviceSynchronize());
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
  no_zero_size = thrust::count_if(d_n_circulations_all_loc_ml.begin(), d_n_circulations_all_loc_ml.end(), NotZero<unsigned int>());
  auto loc_index_begin_2 = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.begin(),
                                                                        d_ce_all_location_target_indices.begin(),
                                                                        d_ce_all_moving_levels.begin()));
  auto loc_index_end_2 = thrust::make_zip_iterator(thrust::make_tuple(d_ce_all_location_from_indices.end(),
                                                                      d_ce_all_location_target_indices.end(),
                                                                      d_ce_all_moving_levels.end()));
  ThrustTuple3VectorDevice<int, int, int> d_circulations_all_loc_ml_indices_no_zero(no_zero_size);
  thrust::copy_if(thrust::device,
                  loc_index_begin_2,
                  loc_index_end_2,
                  d_n_circulations_all_loc_ml.begin(),
                  d_circulations_all_loc_ml_indices_no_zero.begin(),
                  CopyNotZero());
  nend = thrust::remove(thrust::device, d_n_circulations_all_loc_ml.begin(), d_n_circulations_all_loc_ml.end(), 0);
  ThrustTVectorDevice<unsigned int> d_n_circulations_all_loc_ml_no_zero(no_zero_size);
  thrust::copy(d_n_circulations_all_loc_ml.begin(), nend, d_n_circulations_all_loc_ml_no_zero.begin());

  /*
   * d_circulations_all_loc_ml_indices_no_zero is tuple of from_location, target_location and moving level
   * d_n_circulations_all_loc_ml_no_zero is non-zero circulation number in all locations and moving levels
   * d_n_circulations_all_loc_ml_no_zero and d_circulations_all_loc_ml_indices_no_zero are the same size and order
   * d_circulate_all_loc_ml_today is tuple of from_location, target_location, moving_level, person_index
   * */
  int total_circulations = thrust::reduce(thrust::device, d_n_circulations_all_loc_ml_no_zero.begin(), d_n_circulations_all_loc_ml_no_zero.end());
  ThrustTuple4VectorDevice<int, int, int, unsigned int> d_circulate_all_loc_ml_today(d_circulations_all_loc_ml_indices_no_zero.size(),
                                                                                     thrust::make_tuple(-1, -1, -1, 0));
  block_size = (d_circulate_all_loc_ml_today.size() + n_threads - 1) / n_threads;
  zip_location_indices_and_n_circulations<<<block_size, n_threads>>>(d_circulate_all_loc_ml_today.size(),
                                                                     thrust::raw_pointer_cast(d_circulations_all_loc_ml_indices_no_zero.data()),
                                                                     thrust::raw_pointer_cast(d_n_circulations_all_loc_ml_no_zero.data()),
                                                                     thrust::raw_pointer_cast(d_circulate_all_loc_ml_today.data()));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());

  thrust::sort(thrust::device, d_circulate_all_loc_ml_today.begin(), d_circulate_all_loc_ml_today.end(), CirculateLessThan());

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
  ThrustTuple5VectorDevice<int, int, int, unsigned int, int> d_circulate_person_indices_today(total_circulations, thrust::make_tuple(-1, -1, -1, -1, -1));
  block_size = (d_n_circulations_all_loc_ml_no_zero.size() + n_threads - 1) / n_threads;
  fill_circulate_person_indices<<<block_size, n_threads>>>(d_circulate_all_loc_ml_today.size(),
                                                           Model::CONFIG->circulation_info().number_of_moving_levels,
                                                           Model::GPU_RANDOM->d_states,
                                                           thrust::raw_pointer_cast(d_circulate_all_loc_ml_today.data()),
                                                           thrust::raw_pointer_cast(d_ce_popsize_by_moving_level.data()),
                                                           thrust::raw_pointer_cast(d_circulate_person_indices_today.data()));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());

  /*
   * Check if there is any index with n_persons > 1 and random on CPU
   * otherwise schedule person events
   * */
  ThrustTuple5VectorHost<int, int, int, unsigned int, int> h_circulate_person_indices_today = d_circulate_person_indices_today;
  auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationMovingLevel>();
  for (int i = 0; i < d_circulate_all_loc_ml_today.size(); i++) {
    int from_location = h_circulate_person_indices_today[i].get<0>();
    int target_location = h_circulate_person_indices_today[i].get<1>();
    int moving_level = h_circulate_person_indices_today[i].get<2>();
    int n_persons = h_circulate_person_indices_today[i].get<3>();
    auto size = static_cast<int>(pi->vPerson()[from_location][moving_level].size());
    if (size == 0) continue;
    if (n_persons == 1) {
      int p_index = h_circulate_person_indices_today[i].get<4>();
      GPU::Person *p = pi->vPerson()[from_location][moving_level][p_index];
      assert(p->host_state() != GPU::Person::DEAD);
      p->today_target_locations()->push_back(target_location);
      p->randomly_choose_target_location();
//            printf("i %d GPU from %d to %d moving level %d n_persons %d p_index %d\n",
//                   i,
//                   from_location,target_location,moving_level,
//                   n_persons,p_index);
    } else {
      for (int j = 0; j < n_persons; j++) {
        int p_index = Model::RANDOM->random_uniform(size);
        GPU::Person *p = pi->vPerson()[from_location][moving_level][p_index];
        assert(p->host_state() != GPU::Person::DEAD);
        p->today_target_locations()->push_back(target_location);
        p->randomly_choose_target_location();
//                printf("i %d j %d CPU from %d to %d moving level %d n_persons %d p_index %d\n",
//                       i,j,
//                       from_location,target_location,moving_level,
//                       n_persons,p_index);
      }
    }
  }

  if (Model::CONFIG->debug_config().enable_debug_text) {
    auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[GPU Population] Update population circulation GPU (" << d_circulations_indices_no_zero.size() << " " << d_num_leavers_from_target_no_zero.size()
      << " " << total_leavers << " " << total_circulations << ") event time: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms ";
  }
}

/* Update relative biting rate
 * This is equal to GPU::Person::get_age_dependent_biting_factor()
 * */
__device__ double get_age_dependent_biting_factor(int age, int birthday, int current_time, int days_in_year) {
  //
  // 0.00 - 0.25  -  6.5
  // 0.25 - 0.50  -  8.0
  // 0.50 - 0.75  -  9.0
  // 0.75 - 1.00  -  9.5
  // 1.00 - 2.00  -  11.0
  // 2.00 - 3.00  -  13.5
  // 3.00 - 4.00  -  15.5
  // 4.00 - 5.00  -  17.5
  // + 2.75kg until 20
  // then divide by 61.5

  if (age < 1) {
    const auto age = ((current_time - birthday) % days_in_year)
                     / static_cast<double>(days_in_year);
    if (age < 0.25) return 0.106;
    if (age < 0.5) return 0.13;
    if (age < 0.75) return 0.1463;
    return 0.1545;
  }
  if (age < 2) return 0.1789;
  if (age < 3) return 0.2195;
  if (age < 4) return 0.2520;
  if (age < 20) return (17.5 + (age - 4) * 2.75) / 61.5;
  return 1.0;
}

struct AddNewGenotypes : public thrust::unary_function<GPU::PersonUpdateInfo, void> {
    GPU::GenotypeDatabase *gpu_genotype_db_;
    Config *config_;

    __host__
    AddNewGenotypes(GPU::GenotypeDatabase *gpu_genotype_db,
                      Config *config) :
        gpu_genotype_db_(gpu_genotype_db),
        config_(config) {}

    __host__
    void operator()(GPU::PersonUpdateInfo x) {
      if (x.parasites_genotype_mutated_number > 0) {
        for (int p_index = 0; p_index < x.parasites_size; p_index++) {
          if (x.parasite_id[p_index] != -1) {
            gpu_genotype_db_->get_genotype(x.parasite_genotype[p_index], config_);
          }
        }
      }
    }
};

__global__ void update_all_individuals_kernel_stream(int offset,
                                                     int size,
                                                     int current_time,
                                                     curandState *d_state,
                                                     ParasiteDensityLevel h_parasite_density_level,
                                                     ImmuneSystemInformation *d_immune_system_information,
                                                     GPU::DrugType::ResistantAALocation *d_drug_res_aa_loc,
                                                     int *d_drug_res_aa_loc_index,
                                                     double mutation_probability_by_locus,
                                                     char *d_gen_mutation_mask,
                                                     int *d_gen_gene_size,
                                                     int *d_gen_max_copies,
                                                     int *d_gen_aa_size,
                                                     int *d_gen_aa_int,
                                                     int *d_gen_aa_int_start_index,
                                                     GPU::PersonUpdateInfo *d_person_update_info) {
  int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
  if (index < offset + size) {
    curandState local_state = d_state[index];
    if (d_person_update_info[index].person_latest_update_time == current_time) return;
    /* Parasite */
    for (int p_index = 0; p_index < d_person_update_info[index].parasites_size; p_index++) {
      if (d_person_update_info[index].parasite_id[p_index] != -1) {
        int duration = current_time - d_person_update_info[index].person_latest_update_time;
        if (d_person_update_info[index].parasite_update_function_type[p_index] == 1) {
          d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index] = h_parasite_density_level.log_parasite_density_asymptomatic;
        }
        if (d_person_update_info[index].parasite_update_function_type[p_index] == 2) {
          double temp = d_immune_system_information->c_max * (1 - d_person_update_info[index].immune_system_component_latest_value)
                        + d_immune_system_information->c_min * d_person_update_info[index].immune_system_component_latest_value;
          d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index] =
              d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index]
              + duration * (log10(temp) + log10(d_person_update_info[index].parasite_genotype_fitness_multiple_infection[p_index]));
        }
      }
    }
    /* Drug */
    if (d_person_update_info[index].person_has_drug_in_blood) {
      for (int d_index = 0; d_index < d_person_update_info[index].drug_in_blood_size; d_index++) {
        const int d_type_id = d_person_update_info[index].drug_in_blood_type_id[d_index];
        if (d_type_id >= 0 && d_type_id < d_person_update_info[index].drug_in_blood_size) {
          d_person_update_info[index].drug_last_update_time[d_type_id] = current_time;
          const auto days = current_time - d_person_update_info[index].drug_start_time[d_type_id];
          if (days == 0) {
            d_person_update_info[index].drug_last_update_value[d_type_id] = 0;
          } else if (days <= d_person_update_info[index].drug_dosing_days[d_type_id]) {
            if (d_type_id == 0) {
              // drug is artemisinin
              d_person_update_info[index].drug_last_update_value[d_type_id] =
                  d_person_update_info[index].drug_starting_value[d_type_id] + d_person_update_info[index].drug_rand_uniform_1[d_type_id];
            } else {
              d_person_update_info[index].drug_starting_value[d_type_id] += days >= 1 ? d_person_update_info[index].drug_rand_uniform_2[d_type_id] : 0;
              d_person_update_info[index].drug_last_update_value[d_type_id] = d_person_update_info[index].drug_starting_value[d_type_id];
            }
          } else {
            const auto temp = is_equal(d_person_update_info[index].drug_half_life[d_type_id], 0.0, d_person_update_info[index].limit_epsilon)
                              ? -100
                              : -(days - d_person_update_info[index].drug_dosing_days[d_type_id]) * logf(2)
                                / d_person_update_info[index].drug_half_life[d_type_id];  //-ai*t = - t* ln2 / tstar
            if (exp(temp) <= (10.0 / 100.0)) {
              d_person_update_info[index].drug_last_update_value[d_type_id] = 0;
            } else {
              d_person_update_info[index].drug_last_update_value[d_type_id] = d_person_update_info[index].drug_starting_value[d_type_id] * exp(temp);
            }
          }
        }
      }
    }
    /* Parasite by drug */
    d_person_update_info[index].parasites_genotype_mutated_number = 0;
    for (int p_index = 0; p_index < d_person_update_info[index].parasites_size; p_index++) {
      if (d_person_update_info[index].parasite_id[p_index] != -1) {
        char *new_gen_aa = d_person_update_info[index].parasite_genotype[p_index];
        if (d_person_update_info[index].person_has_drug_in_blood) {
          for (int d_index = 0; d_index < d_person_update_info[index].drug_in_blood_size; d_index++) {
            const int d_type_id = d_person_update_info[index].drug_in_blood_type_id[d_index];
            if (d_type_id >= 0 && d_drug_res_aa_loc_index[d_type_id] != -1
                && d_type_id < d_person_update_info[index].drug_in_blood_size) {
              int *d_drug_res_aa_loc_index_int = (int *) malloc(2 * sizeof(int));
              decode_int_to_arr2(d_drug_res_aa_loc_index[d_type_id], d_drug_res_aa_loc_index_int);
              for (int aa_pos_id = d_drug_res_aa_loc_index_int[0];
                   aa_pos_id < d_drug_res_aa_loc_index_int[1];
                   aa_pos_id++) {
                if (d_gen_mutation_mask[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string] == '1') {
                  const double p = curand_gsl_uniform_double(curand_uniform_double(&local_state), 0.0, 1.0);
                  if (p < mutation_probability_by_locus) {
                    if (d_drug_res_aa_loc[aa_pos_id].is_copy_number) {
                      int max_copies = -1;
                      if (d_gen_gene_size[d_drug_res_aa_loc[aa_pos_id].chromosome_id] > 1) {
                        int *max_copies_int = (int *) malloc(2 * sizeof(int));
                        decode_int_to_arr2(d_gen_max_copies[d_drug_res_aa_loc[aa_pos_id].chromosome_id], max_copies_int);
                        max_copies = max_copies_int[d_drug_res_aa_loc[aa_pos_id].gene_id];
                        free(max_copies_int);
                      } else {
                        max_copies = d_gen_max_copies[d_drug_res_aa_loc[aa_pos_id].chromosome_id];
                      }
                      int old_copy_number = (int) (new_gen_aa[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string]) - 48;
                      if (old_copy_number == 1) {
                        new_gen_aa[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string] = char((old_copy_number + 1) + 48);
                      } else if (old_copy_number == max_copies) {
                        new_gen_aa[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string] = char((old_copy_number - 1) + 48);
                      } else {
                        int new_copy_number = curand_gsl_uniform_double(curand_uniform_double(&local_state), 0.0, 1.0) < 0.5 ?
                                              old_copy_number - 1 : old_copy_number + 1;
                        new_gen_aa[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string] = char((new_copy_number) + 48);
                      }
                    } else {
                      int aa_start_index = d_gen_aa_int_start_index[d_drug_res_aa_loc[aa_pos_id].chromosome_id];
                      int aa_int_index = aa_start_index + d_drug_res_aa_loc[aa_pos_id].aa_id;
                      int *aa_list_int = (int *) malloc(2 * sizeof(int));
                      decode_int_to_arr2(d_gen_aa_int[aa_int_index], aa_list_int);
                      int aa_list_size = sizeof(aa_list_int) / sizeof(int);
                      char old_aa = char(aa_list_int[d_drug_res_aa_loc[aa_pos_id].gene_id] + 48);
                      // draw random aa id
                      int new_aa_id = curand_gsl_uniform_int(curand_uniform_double(&local_state), 0, aa_list_size - 1);
                      char new_aa = char(aa_list_int[new_aa_id] + 48);
                      if (new_aa == old_aa) {
                        new_aa = char(aa_list_int[new_aa_id + 1] + 48);
                      }
                      new_gen_aa[d_drug_res_aa_loc[aa_pos_id].aa_index_in_aa_string] = new_aa;
                      free(aa_list_int);
                    }
                    d_person_update_info[index].parasites_genotype_mutated_number += 1;
                  }
                }
              }
              free(d_drug_res_aa_loc_index_int);
            }
          }
        }
      }
    }
    /* Immune */
    auto immune_component_temp = 0.0;
    if (d_person_update_info[index].immune_system_is_increased) {
      //increase I(t) = 1 - (1-I0)e^(-b1*t)
      double immune_component_acquire_rate = 0.0;
      if (d_person_update_info[index].immune_system_component_type == 1) {
        /* from InfantImmuneComponent acquire */
        immune_component_acquire_rate = 0.0;
      }
      if (d_person_update_info[index].immune_system_component_type == 2) {
        /* from NonInfantImmuneComponent acquire */
        immune_component_acquire_rate = (d_person_update_info[index].person_age > 80)
                                        ? d_immune_system_information->acquire_rate_by_age[80]
                                        : d_immune_system_information->acquire_rate_by_age[d_person_update_info[index].person_age];
      }
      immune_component_temp = 1.0 - (1.0 - d_person_update_info[index].immune_system_component_latest_value)
                                    * exp(-immune_component_acquire_rate * (current_time - d_person_update_info[index].person_latest_update_time));
    } else {
      //decrease I(t) = I0 * e ^ (-b2*t);
      double immune_component_decay_rate = 0.0;
      if (d_person_update_info[index].immune_system_component_type == 1) {
        /* from InfantImmuneComponent decay */
        immune_component_decay_rate = 0.0315;
      }
      if (d_person_update_info[index].immune_system_component_type == 2) {
        /* from NonInfantImmuneComponent decay */
        immune_component_decay_rate = d_immune_system_information->decay_rate;
      }
      immune_component_temp = d_person_update_info[index].immune_system_component_latest_value
                              * exp(-immune_component_decay_rate * (current_time - d_person_update_info[index].person_latest_update_time));
      immune_component_temp = (immune_component_temp < 0.00001) ? 0.0 : immune_component_temp;
    }
    d_person_update_info[index].immune_system_component_latest_value = immune_component_temp;
    /* Cutoff */
    if (d_person_update_info[index].person_has_drug_in_blood) {
      for (int d_index = 0; d_index < d_person_update_info[index].drug_in_blood_size; d_index++) {
        const int d_type_id = d_person_update_info[index].drug_in_blood_type_id[d_index];
        if (d_type_id >= 0 && d_type_id < d_person_update_info[index].drug_in_blood_size) {
          if (d_person_update_info[index].drug_last_update_value[d_type_id] <= DRUG_CUT_OFF_VALUE) {
            d_person_update_info[index].drug_in_blood_type_id[d_type_id] = -1;
          }
        }
      }
    }
    /* Clear cured */
    for (int p_index = 0; p_index < d_person_update_info[index].parasites_size; p_index++) {
      if (d_person_update_info[index].parasite_id[p_index] != -1) {
        if (d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index]
            <= h_parasite_density_level.log_parasite_density_cured + 0.00001) {
          d_person_update_info[index].parasite_id[p_index] = -1;
          d_person_update_info[index].parasite_update_function_type[p_index] = 0;
          d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index] = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
          d_person_update_info[index].parasite_genotype_fitness_multiple_infection[p_index] = 1.0;
          d_person_update_info[index].parasite_gametocyte_level[p_index] = 0.0;
          d_person_update_info[index].parasite_log10_infectious_density[p_index] = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
          d_person_update_info[index].parasites_current_index -= 1;
          d_person_update_info[index].parasites_size -= 1;
        } else {
          /* From GPU::ClonalParasitePopulation::get_log10_infectious_density() */
          if (is_equal(d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index],
                       d_person_update_info[index].LOG_ZERO_PARASITE_DENSITY,
                       d_person_update_info[index].limit_epsilon)
              || is_equal(d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index],
                          0.0,
                          d_person_update_info[index].limit_epsilon)) {
            d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index]
                = d_person_update_info[index].LOG_ZERO_PARASITE_DENSITY;
          }
          d_person_update_info[index].parasite_log10_infectious_density[p_index]
              = d_person_update_info[index].parasite_last_update_log10_parasite_density[p_index]
                + log10(d_person_update_info[index].parasite_gametocyte_level[p_index]);

          /* From GPU::SingleHostClonalParasitePopulations::clear_cured_parasites() */
          if (d_person_update_info[index].parasites_log10_total_infectious_density == d_person_update_info[index].LOG_ZERO_PARASITE_DENSITY) {
            d_person_update_info[index].parasites_log10_total_infectious_density
                = d_person_update_info[index].parasite_log10_infectious_density[p_index];
          } else {
            d_person_update_info[index].parasites_log10_total_infectious_density
                += log10(pow(10, d_person_update_info[index].parasite_log10_infectious_density[p_index]
                                 - d_person_update_info[index].parasites_log10_total_infectious_density) + 1);
          }
        }
      }
    }
    /* Change state */
    if (d_person_update_info[index].parasites_size == 0) {
      if (d_person_update_info[index].person_liver_parasite_genotype[0] == '\0') {
        d_person_update_info[index].person_host_state = static_cast<int>(GPU::Person::SUSCEPTIBLE);
      } else {
        d_person_update_info[index].person_host_state = static_cast<int>(GPU::Person::EXPOSED);
      }
      d_person_update_info[index].immune_system_is_increased = false;
    } else {
      d_person_update_info[index].immune_system_is_increased = true;
    }
    if (d_person_update_info[index].person_using_age_dependent_biting_level) {
      d_person_update_info[index].person_current_relative_biting_rate
          = d_person_update_info[index].person_innate_relative_biting_rate
            * get_age_dependent_biting_factor(d_person_update_info[index].person_age,
                                              d_person_update_info[index].person_birthday,
                                              current_time,
                                              365);
    } else {
      d_person_update_info[index].person_current_relative_biting_rate
          = d_person_update_info[index].person_innate_relative_biting_rate;
    }
    /* Latest update time */
    d_person_update_info[index].person_latest_update_time = current_time;
//        if(index == offset + size - 1) printf("GPU all_individuals_kernel_stream (%d -> %d) index %d last_time %d curr_time %d\n",
//                                              offset,offset+size,index,d_person_update_info[index].person_latest_update_time,current_time);
    d_state[index] = local_state;
    __syncthreads();
  }
}

void GPU::PopulationKernel::update_all_individuals() {
//  for(int i = pi->h_person_update_info().size() - 5; i < pi->h_person_update_info().size(); i++){
//    LOG_IF(Model::CONFIG->debug_config().enable_debug_text,INFO)
//      << fmt::format("[PopulationKernel update_all_individuals] BEFORE Person {} last update time {} parasites_size {} parasites_genotype_mutated_number {}"
//                     " person_biting {} person_moving {}",
//                     i,
//                     pi->h_person_update_info()[i].person_latest_update_time,
//                     pi->h_person_update_info()[i].parasites_size,
//                     pi->h_person_update_info()[i].parasites_genotype_mutated_number,
//                     pi->h_person_update_info()[i].person_current_relative_biting_rate,
//                     pi->h_person_update_info()[i].person_current_relative_moving_rate);
//  }
  check_cuda_error(cudaEventRecord(start_event, 0));
  LOG_IF(Model::CONFIG->debug_config().enable_debug_text
         && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
    << fmt::format("[PopulationKernel update_all_individuals] Start size {} {}GBs",
                   pi->h_person_update_info().size(),
                   pi->h_person_update_info().size() * sizeof(GPU::PersonUpdateInfo) / 1e9);
  int batch_size = (pi->h_person_update_info().size() < Model::CONFIG->gpu_config().n_people_1_batch)
                   ? pi->h_person_update_info().size() : Model::CONFIG->gpu_config().n_people_1_batch;
  int batch_offset = 0;
  for (auto [batch_remain, b_index] = std::tuple{pi->h_person_update_info().size(), 0}; batch_remain > 0; batch_remain -= batch_size, b_index++) {
    batch_size = (batch_remain < batch_size) ? batch_remain : batch_size;
    int batch_from = pi->h_person_update_info().size() - batch_remain;
    int batch_to = batch_from + batch_size;
    const int batch_bytes = batch_size * sizeof(GPU::PersonUpdateInfo);
    check_cuda_error(cudaMalloc((void **) &d_buffer_person_update_info_stream, batch_bytes));
//    LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//           && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//           << fmt::format("[PopulationKernel update_all_individuals] Work batch size {} remain {}, from {} to {} (of {})",
//                                batch_size, batch_remain, batch_from, batch_to,
//                                pi->h_person_update_info().size());
    /*
     * Check if batch_size > Model::CONFIG->gpu_config().n_people_1_batch / 2
     * If true, then use streams
     * */
    if(batch_size > (Model::CONFIG->gpu_config().n_people_1_batch / 2)){
      int stream_size = batch_size / Model::CONFIG->gpu_config().n_streams;
      int stream_offset = 0;
      for (auto [stream_remain, s_index] = std::tuple{batch_size, 0}; stream_remain > 0; stream_remain -= stream_size, s_index++) {
        if (s_index == Model::CONFIG->gpu_config().n_streams - 1) {
          stream_size = stream_remain;
        } else {
          stream_size = (stream_remain < stream_size) ? stream_remain : stream_size;
        }
        const int stream_bytes = stream_size * sizeof(GPU::PersonUpdateInfo);
        const int n_threads = (stream_size < Model::CONFIG->gpu_config().n_threads) ? stream_size : Model::CONFIG->gpu_config().n_threads;
        const int n_blocks = (stream_size + n_threads - 1) / n_threads;
//        LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//               && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//          << fmt::format("[PopulationKernel update_all_individuals]   Work stream {} size {}, from {} to {} offset {} n_threads {} n_blocks {} per stream",
//                         s_index, stream_size, batch_size - stream_remain, batch_size - stream_remain + stream_size,
//                         stream_offset,
//                         n_threads,
//                         n_blocks);
        check_cuda_error(cudaMemcpyAsync(&d_buffer_person_update_info_stream[stream_offset],
                                         thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                         stream_bytes,
                                         cudaMemcpyHostToDevice,
                                         d_streams[s_index]));
        update_all_individuals_kernel_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(stream_offset,
                                                                                             stream_size,
                                                                                             Model::GPU_SCHEDULER->current_time(),
                                                                                             Model::GPU_RANDOM->d_states,
                                                                                             Model::CONFIG->parasite_density_level(),
                                                                                             d_immune_system_information,
                                                                                             thrust::raw_pointer_cast(d_drug_res_aa_loc.data()),
                                                                                             thrust::raw_pointer_cast(d_drug_res_aa_loc_index.data()),
                                                                                             Model::CONFIG->mutation_probability_by_locus(),
                                                                                             d_gen_mutation_mask,
                                                                                             thrust::raw_pointer_cast(d_gen_gene_size.data()),
                                                                                             thrust::raw_pointer_cast(d_gen_max_copies.data()),
                                                                                             thrust::raw_pointer_cast(d_gen_aa_size.data()),
                                                                                             thrust::raw_pointer_cast(d_gen_aa_int.data()),
                                                                                             thrust::raw_pointer_cast(d_gen_aa_int_start_index.data()),
                                                                                             d_buffer_person_update_info_stream);
        check_cuda_error(cudaDeviceSynchronize());
        check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                         &d_buffer_person_update_info_stream[stream_offset],
                                         stream_bytes,
                                         cudaMemcpyDeviceToHost,
                                         d_streams[s_index]));
        check_cuda_error(cudaGetLastError());
        stream_offset += stream_size;
        batch_offset += stream_size;
      }
      for (int s_index = 0; s_index < Model::CONFIG->gpu_config().n_streams; s_index++) {
        check_cuda_error(cudaStreamSynchronize(d_streams[s_index]));
      }
    }
    else{
      const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
      const int n_blocks = (batch_size + n_threads - 1) / n_threads;
      int s_index = 0;
//      LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//             && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//        << fmt::format("[PopulationKernel update_all_individuals]   Last batch {} size {}, from {} to {} offset {} n_threads {} n_blocks {} per batch",
//                       s_index, batch_size, batch_from, batch_to,
//                       batch_offset,
//                       n_threads,
//                       n_blocks);
      check_cuda_error(cudaMemcpyAsync(&d_buffer_person_update_info_stream[0],
                                       thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                       batch_bytes,
                                       cudaMemcpyHostToDevice,
                                       d_streams[s_index]));
      update_all_individuals_kernel_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(0,
                                                                                           batch_size,
                                                                                           Model::GPU_SCHEDULER->current_time(),
                                                                                           Model::GPU_RANDOM->d_states,
                                                                                           Model::CONFIG->parasite_density_level(),
                                                                                           d_immune_system_information,
                                                                                           thrust::raw_pointer_cast(d_drug_res_aa_loc.data()),
                                                                                           thrust::raw_pointer_cast(d_drug_res_aa_loc_index.data()),
                                                                                           Model::CONFIG->mutation_probability_by_locus(),
                                                                                           d_gen_mutation_mask,
                                                                                           thrust::raw_pointer_cast(d_gen_gene_size.data()),
                                                                                           thrust::raw_pointer_cast(d_gen_max_copies.data()),
                                                                                           thrust::raw_pointer_cast(d_gen_aa_size.data()),
                                                                                           thrust::raw_pointer_cast(d_gen_aa_int.data()),
                                                                                           thrust::raw_pointer_cast(d_gen_aa_int_start_index.data()),
                                                                                           d_buffer_person_update_info_stream);
      check_cuda_error(cudaDeviceSynchronize());
      check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                       &d_buffer_person_update_info_stream[0],
                                       batch_bytes,
                                       cudaMemcpyDeviceToHost,
                                       d_streams[s_index]));
      check_cuda_error(cudaGetLastError());
      batch_offset += batch_size;
    }
    check_cuda_error(cudaFree(d_buffer_person_update_info_stream));
    /*
     * Add new genotypes to genotype database on host
     * First copy new genotype to new device vector,
     * then copy to host
     * */
    thrust::for_each(thrust::host,
                     pi->h_person_update_info().begin() + batch_from,
                     pi->h_person_update_info().begin() + batch_to,
                     AddNewGenotypes(&Model::CONFIG->gpu_genotype_db,
                                       Model::CONFIG));
    LOG_IF(Model::CONFIG->debug_config().enable_debug_text
           && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
      << fmt::format("[PopulationKernel update_all_individuals] Batch size {} gpu_genotype_db size {}",
                     batch_size, Model::CONFIG->gpu_genotype_db.size());
  }
  check_cuda_error(cudaEventRecord (stop_event, 0));
  check_cuda_error(cudaEventSynchronize (stop_event));
  check_cuda_error(cudaEventElapsedTime (&elapsed_time, start_event, stop_event));
  LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
  << fmt::format("[PopulationKernel::update_all_individuals] Update population all individuals GPU time: {} ms", elapsed_time);
//  for(int i = pi->h_person_update_info().size() - 5; i < pi->h_person_update_info().size(); i++){
//    LOG_IF(Model::CONFIG->debug_config().enable_debug_text,INFO)
//      << fmt::format("[PopulationKernel update_all_individuals] AFTER Person {} last update time {} parasites_size {} parasites_genotype_mutated_number {}"
//                     " person_biting {} person_moving {}",
//                     i,
//                     pi->h_person_update_info()[i].person_latest_update_time,
//                     pi->h_person_update_info()[i].parasites_size,
//                     pi->h_person_update_info()[i].parasites_genotype_mutated_number,
//                     pi->h_person_update_info()[i].person_current_relative_biting_rate,
//                     pi->h_person_update_info()[i].person_current_relative_moving_rate);
//  }
}

__global__ void update_current_foi_kernel_stream(int offset,
                                                 int size,
                                                 RelativeInfectivity h_relative_infectivity,
                                                 GPU::PersonUpdateInfo *d_person_update_info) {
  int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
  if (index < offset + size) {
    /* relative infectivity and FOI */
    if (d_person_update_info[index].parasites_log10_total_infectious_density == d_person_update_info[index].LOG_ZERO_PARASITE_DENSITY ||
        d_person_update_info[index].person_host_state == 4) {
      d_person_update_info[index].person_current_relative_infectivity = 0.0;
    }
    else{
      // this sigma has already taken 'ln' and 'log10' into account
      const auto d_n = d_person_update_info[index].parasites_log10_total_infectious_density  * h_relative_infectivity.sigma
                       + h_relative_infectivity.ro_star;
      const auto p = 0.001;//curand_gsl_cdf_ugaussian_P(d_n);
      d_person_update_info[index].person_current_relative_infectivity = p * p + 0.01;
    }
    d_person_update_info[index].person_current_foi = d_person_update_info[index].person_current_relative_biting_rate
                                                     * d_person_update_info[index].person_current_relative_infectivity;
    __syncthreads();
  }
}

struct LocationLessThan {
    __host__ __device__
    bool operator()(thrust::tuple<int, double, double, double> t0, thrust::tuple<int, double, double, double> t1) {
      if (thrust::get<0>(t0) < thrust::get<0>(t1))
        return true;
      if (thrust::get<0>(t0) > thrust::get<0>(t1))
        return false;
      return thrust::get<0>(t0) < thrust::get<0>(t1);
    }

};

struct CheckKeyTuple {
    __host__ __device__
    bool operator()(ThrustTuple4<int, double, double, double> t0, ThrustTuple4<int, double, double, double> t1) {
      if (thrust::get<0>(t0) == thrust::get<0>(t1)) return true;
      else return false;
    }
};

struct SumValueTuple {
    __host__ __device__
    ThrustTuple4<int, double, double, double> operator()(ThrustTuple4<int, double, double, double> t0, ThrustTuple4<int, double, double, double> t1) {
      if(thrust::get<0>(t0) == thrust::get<0>(t1)){
        return thrust::make_tuple(thrust::get<0>(t0), thrust::get<1>(t0) + thrust::get<1>(t1),
                                  thrust::get<2>(t0) + thrust::get<2>(t1), thrust::get<3>(t0) + thrust::get<3>(t1));
      }
    }
};

void GPU::PopulationKernel::update_current_foi() {
  check_cuda_error(cudaEventRecord(start_event, 0));
  d_sum_biting_moving_foi_by_loc.resize(0);
  int batch_size = (pi->h_person_update_info().size() < Model::CONFIG->gpu_config().n_people_1_batch)
                   ? pi->h_person_update_info().size() : Model::CONFIG->gpu_config().n_people_1_batch;
  int batch_offset = 0;
  for (auto [batch_remain, b_index] = std::tuple{pi->h_person_update_info().size(), 0}; batch_remain > 0; batch_remain -= batch_size, b_index++) {
    batch_size = (batch_remain < batch_size) ? batch_remain : batch_size;
    int batch_from = pi->h_person_update_info().size() - batch_remain;
    int batch_to = batch_from + batch_size;
    const int batch_bytes = batch_size * sizeof(GPU::PersonUpdateInfo);
    check_cuda_error(cudaMalloc((void **) &d_buffer_person_update_info_stream, batch_bytes));
//    LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//           && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//      << fmt::format("[PopulationKernel update_current_foi] Work batch size {} remain {}, from {} to {} (of {})",
//                     batch_size, batch_remain, batch_from, batch_to,
//                     pi->h_person_update_info().size());
    /*
     * Check if batch_size > Model::CONFIG->gpu_config().n_people_1_batch / 2
     * If true, then use streams
     * */
    if(batch_size > (Model::CONFIG->gpu_config().n_people_1_batch / 2)){
      int stream_size = batch_size / Model::CONFIG->gpu_config().n_streams;
      int stream_offset = 0;
      for (auto [stream_remain, s_index] = std::tuple{batch_size, 0}; stream_remain > 0; stream_remain -= stream_size, s_index++) {
        if (s_index == Model::CONFIG->gpu_config().n_streams - 1) {
          stream_size = stream_remain;
        } else {
          stream_size = (stream_remain < stream_size) ? stream_remain : stream_size;
        }
        const int stream_bytes = stream_size * sizeof(GPU::PersonUpdateInfo);
        const int n_threads = (stream_size < Model::CONFIG->gpu_config().n_threads) ? stream_size : Model::CONFIG->gpu_config().n_threads;
        const int n_blocks = (stream_size + n_threads - 1) / n_threads;
//        LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//               && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//          << fmt::format("[PopulationKernel update_current_foi]   Work stream {} size {}, from {} to {} b_offset {} s_offset {} n_threads {} n_blocks {} per stream",
//                         s_index, stream_size, batch_size - stream_remain, batch_size - stream_remain + stream_size,
//                         batch_offset,
//                         stream_offset,
//                         n_threads,
//                         n_blocks);
        check_cuda_error(cudaMemcpyAsync(&d_buffer_person_update_info_stream[stream_offset],
                                         thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                         stream_bytes,
                                         cudaMemcpyHostToDevice,
                                         d_streams[s_index]));
        update_current_foi_kernel_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(stream_offset,
                                                                                         stream_size,
                                                                                         Model::CONFIG->relative_infectivity(),
                                                                                         d_buffer_person_update_info_stream
                                                                                         );
        check_cuda_error(cudaDeviceSynchronize());
        check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                         &d_buffer_person_update_info_stream[stream_offset],
                                         stream_bytes,
                                         cudaMemcpyDeviceToHost,
                                         d_streams[s_index]));
        check_cuda_error(cudaGetLastError());
        stream_offset += stream_size;
        batch_offset += stream_size;
      }
      for (int s_index = 0; s_index < Model::CONFIG->gpu_config().n_streams; s_index++) {
        check_cuda_error(cudaStreamSynchronize(d_streams[s_index]));
      }
    }
    else{
      const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
      const int n_blocks = (batch_size + n_threads - 1) / n_threads;
      int s_index = 0;
//      LOG_IF(Model::CONFIG->debug_config().enable_debug_text
//             && Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0,INFO)
//        << fmt::format("[PopulationKernel update_current_foi]   Last batch {} size {}, from {} to {} b_offset {} n_threads {} n_blocks {} per batch",
//                       s_index, batch_size, batch_from, batch_to,
//                       batch_offset,
//                       n_threads,
//                       n_blocks);
      check_cuda_error(cudaMemcpyAsync(&d_buffer_person_update_info_stream[0],
                                       thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                       batch_bytes,
                                       cudaMemcpyHostToDevice,
                                       d_streams[s_index]));
      update_current_foi_kernel_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(0,
                                                                                batch_size,
                                                                                Model::CONFIG->relative_infectivity(),
                                                                                d_buffer_person_update_info_stream);
      check_cuda_error(cudaDeviceSynchronize());
      check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_update_info().data()) + batch_offset,
                                       &d_buffer_person_update_info_stream[0],
                                       batch_bytes,
                                       cudaMemcpyDeviceToHost,
                                       d_streams[s_index]));
      check_cuda_error(cudaGetLastError());
      batch_offset += batch_size;
    }
    /*
     * Sum up biting, moving and foi by location after every batch
     * The length of d_sum_biting_moving_foi_by_loc is 2*number_of_locations
     * */
    auto result = Model::GPU_UTILS->device_sum_biting_moving_foi_by_loc_pointer(d_buffer_person_update_info_stream,0,batch_size);
    d_sum_biting_moving_foi_by_loc.reserve(Model::CONFIG->number_of_locations()*2);
    d_sum_biting_moving_foi_by_loc.insert(d_sum_biting_moving_foi_by_loc.end(),result.begin(),result.begin()+Model::CONFIG->number_of_locations());
    h_sum_biting_moving_foi_by_loc = Model::GPU_UTILS->host_reduce_vector(d_sum_biting_moving_foi_by_loc);
    check_cuda_error(cudaDeviceSynchronize());
    check_cuda_error(cudaGetLastError());
    check_cuda_error(cudaFree(d_buffer_person_update_info_stream));
  }
  check_cuda_error(cudaEventRecord (stop_event, 0));
  check_cuda_error(cudaEventSynchronize (stop_event));
  check_cuda_error(cudaEventElapsedTime (&elapsed_time, start_event, stop_event));
  LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
    << fmt::format("[PopulationKernel::update_current_foi] Update current FOI GPU time: {} ms", elapsed_time);
//  if(Model::CONFIG->debug_config().enable_debug_text) {
//    for (int loc = Model::CONFIG->number_of_locations() - 5; loc < Model::CONFIG->number_of_locations(); loc++) {
//      printf("%d GPU sum_relative_biting_by_location[%d] biting %f\n", Model::GPU_SCHEDULER->current_time(),
//             h_sum_biting_moving_foi_by_loc[loc].get<0>(), h_sum_biting_moving_foi_by_loc[loc].get<1>());
//      printf("%d GPU sum_relative_moving_by_location[%d] moving %f\n", Model::GPU_SCHEDULER->current_time(),
//             h_sum_biting_moving_foi_by_loc[loc].get<0>(), h_sum_biting_moving_foi_by_loc[loc].get<2>());
//      printf("%d GPU sum_relative_moving_by_location[%d] foi %f\n", Model::GPU_SCHEDULER->current_time(),
//             h_sum_biting_moving_foi_by_loc[loc].get<0>(), h_sum_biting_moving_foi_by_loc[loc].get<3>());
//    }
//  }
}

void GPU::PopulationKernel::persist_current_force_of_infection_to_use_N_days_later() {
  auto start = std::chrono::high_resolution_clock::now();
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    force_of_infection_for_N_days_by_location[Model::GPU_SCHEDULER->current_time()
                                              % Model::CONFIG->number_of_tracking_days()][loc] =
    h_sum_biting_moving_foi_by_loc[loc].get<3>();
  }
  auto lapse = std::chrono::high_resolution_clock::now() - start;
  if(Model::CONFIG->debug_config().enable_debug_text){
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[PopulationKernel::persist_current_force_of_infection_to_use_N_days_later] Update FOI N days CPU time: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
}

//
//void GPU::PopulationKernel::calculate_n_person_bitten_today(int n_locations,
//                                                            ThrustTVectorDevice<double> &d_foi_all_locations,
//                                                            ThrustTVectorDevice<int> &d_n_person_bitten_today_all_locations) {
//
//
//}
//
//void GPU::PopulationKernel::perform_infection_event() {
//  auto tp_start = std::chrono::high_resolution_clock::now();
//  auto tracking_index = Model::GPU_SCHEDULER->current_time() % Model::CONFIG->number_of_tracking_days();
//
//  /*
//   * Calculate probability of leaving location in all locations (n_location*n_location)
//   * Also get indices from and to arrays
//   * */
//
//  ThrustTVectorDevice<double> d_foi_all_locations;
//  ThrustTVectorDevice<int> d_n_person_bitten_today_all_locations;
//  calculate_n_person_bitten_today(Model::CONFIG->number_of_locations(),
//                                  d_foi_all_locations, d_n_person_bitten_today_all_locations);
//}



