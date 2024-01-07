//
// Created by kient on 6/17/2023.
//

#include <thrust/shuffle.h>
#include <thrust/random.h>
#include <thrust/execution_policy.h>
#include <thrust/count.h>
#include <thrust/sort.h>
#include "Random.cuh"
#include "Utils.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

GPU::Random::Random() {
    d_states = nullptr;
    n_threads = 1024;
}

GPU::Random::~Random() {
    cudaFree(d_states);
}

__global__ void setup(int num,curandState *state, long seed)
{
    auto id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < num) curand_init(seed, id, 0, &state[id]);
}

struct setupCurandStates
{
    long seed_;
    setupCurandStates(long seed) : seed_(seed) {}
    __device__
    curandState operator()(int id){
        curandState s;
        curand_init(seed_, id, 0, &s);
        return s;
    }
};

void GPU::Random::init_curand_states(ThrustTVectorDevice<curandState> &d_curand_states, int size, long seed){
    thrust::transform(thrust::counting_iterator<int>(0),
                      thrust::counting_iterator<int>(size),
                      d_curand_states.begin(),
                      setupCurandStates(seed));
}

void GPU::Random::init(int n, unsigned long seed) {
    LOG(INFO) << "GPU Random initializing with seed: " << seed;
    cudaMalloc((void **) &d_states, sizeof(curandState) * n);
    n_blocks = (n + n_threads + 1) / n_threads;
    setup<<<n_blocks,n_threads>>>(n,d_states, seed);
    check_cuda_error(cudaDeviceSynchronize());
    check_cuda_error(cudaPeekAtLastError());
}

void GPU::Random::free() {
    cudaFree(d_states);
}

/*
 * Binomial are from
 * https://stackoverflow.com/questions/23561551/a-efficient-binomial-random-number-generator-code-in-java/23574723#23574723
 * to get O(Np) when p is small instead of O(N) in curand_binomial_naive
 * Not need to use Rejection because N is not very large, as mentioned in
 * https://peterchng.com/blog/2020/10/23/building-binomial-and-multinomial-samplers-in-java/
 */
// Function to generate random numbers from a binomial distribution
__device__ unsigned int curand_binomial(curandState *state, double p, unsigned int N)
{
    double log_q = log(1.0 - p);
    int count = 0;
    double sum = 0;
    while(true) {
        sum += log(curand_uniform_double(state)) / (N - count);
        if(sum < log_q) {
            return count;
        }
        count++;
    }
}

/*
 * Binomial naive from
 * https://peterchng.com/blog/2020/10/23/building-binomial-and-multinomial-samplers-in-java/
 */
__device__ unsigned int curand_binomial_naive(curandState *state, double p, unsigned int N)
{
    unsigned int count = 0;
    if(N <= 0) return 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        double rand_num = curand_uniform_double(state);
        if (rand_num < p)
        {
            count++;
        }
    }
    return count;
}

/*
 * Multinomial are from gsl_multinomial
 * This is parallel multinomial on n_locations, with K draws at each location.
 * The output must have size n_locations*K
 * N is the number of trials, p is the probability of each outcome, n is the output
 */
__global__ void multinomial_kernel(curandState *d_state, int n_locations, int K, int N[], double p[], unsigned int n[],
                                   double norm[], double sum_p[], unsigned int sum_n[])
{
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    curandState local_state = d_state[thread_index];
    for(int index = thread_index; index < n_locations; index += stride){
        if(N[index] <= 0) {
            for (size_t i = 0; i < K; i++){
                n[index * K + i] = 0;
            }
            return;
        }
        for (size_t i = 0; i < K; i++)
        {
            norm[index] += p[index*K+i];
        }

        __syncthreads();

        // Calculate multinomial distribution
        for (size_t i = 0; i < K; i++)
        {
            if (p[index*K+i] > 0.0)
            {
                n[index*K+i] = curand_binomial_naive(&local_state, p[index*K+i] / (norm[index] - sum_p[index]), N[index] - sum_n[index]);
            }
            else
            {
                n[index*K+i] = 0;
            }

            sum_p[index] += p[index*K+i];
            sum_n[index] += n[index*K+i];
        }
        __syncthreads();
    }
    d_state[thread_index] = local_state;
}

/*
 * Multinomial are from gsl_multinomial
 * This is parallel multinomial on n_locations, with d_n_trials at each location.
 * The output must have size n_locations*n_samples_each_location
 */
void GPU::Random::random_multinomial(int n_locations, int n_samples_each_location,
                                     ThrustTVectorDevice<int> d_n_trials,
                                     ThrustTVectorDevice<double> d_distributions,
                                     ThrustTVectorDevice<unsigned int> &d_samples){
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int n_blocks = (n_locations + n_threads + 1) / n_threads;
    ThrustTVectorDevice<double> d_norm(n_locations,0.0);
    ThrustTVectorDevice<double> d_sum_p(n_locations,0.0);
    ThrustTVectorDevice<unsigned int> d_sum_n(n_locations,0);
    multinomial_kernel<<<n_blocks, n_threads>>>(d_states,
                                                n_locations,
                                                n_samples_each_location,
                                                thrust::raw_pointer_cast(d_n_trials.data()),
                                                thrust::raw_pointer_cast(d_distributions.data()),
                                                thrust::raw_pointer_cast(d_samples.data()),
                                                thrust::raw_pointer_cast(d_norm.data()),
                                                thrust::raw_pointer_cast(d_sum_p.data()),
                                                thrust::raw_pointer_cast(d_sum_n.data()));
    check_cuda_error(cudaDeviceSynchronize());
    check_cuda_error(cudaPeekAtLastError());
    d_norm.clear();
    d_sum_p.clear();
    d_sum_n.clear();
    ThrustTVectorDevice<double>().swap(d_norm);
    ThrustTVectorDevice<double>().swap(d_sum_p);
    ThrustTVectorDevice<unsigned int>().swap(d_sum_n);
}

/*
 * curand_uniform return (0,1]
 * gsl_rng_uniform return [0,1)
 * */
__global__ void random_uniform_kernel(curandState *d_state,int n_locations, int n_samples_each_location,
                                      double *d_sum_distribution,double *d_uniform_sampling){
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    curandState local_state = d_state[thread_index];
    for(int index = thread_index; index < n_locations*n_samples_each_location; index += stride){
        int location_index = index / n_samples_each_location;
//        printf("kernel index %d location_index %d, sum %f\n",
//               index,location_index,d_sum_distribution[location_index]);
        d_uniform_sampling[index] = curand_uniform_double(&local_state) * d_sum_distribution[location_index];
    }
    d_state[thread_index] = local_state;
}

/*
 * d_n_samples size is n_locations
 * d_distribution_all_locations size is n_locations*n_distributions_each_location
 * d_sum_weight size is n_locations
 * d_all_objects_index size is n_locations*n_samples_each_location
 * d_uniform_sampling size is n_locations*n_samples_each_location
 * d_uniform_sampling_index size is n_locations
 * d_sample_index size is n_locations*n_samples_each_location
 * */
__global__ void roulette_sampling_kernel(int n_locations,
                                         int n_samples_each_location,
                                         int n_distributions_each_location,
                                         double *d_distribution_all_locations,
                                         double *d_sum_weight,
                                         int *d_all_objects_index,
                                         double *d_uniform_sampling,
                                         int *d_uniform_sampling_index,
                                         int *d_sample_index){
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations; index += stride){
        for (auto i = 0; i < n_distributions_each_location; i++) {
            int pi = index * n_distributions_each_location + i;
            d_sum_weight[index] += d_distribution_all_locations[pi];
//            printf("kernel location %d, n_distributions_each_location %d, pi %d, d_distribution_all_locations %f d_uniform_sampling[%d] = %f sum_weight %f\n",
//                   index, n_distributions_each_location, pi, d_distribution_all_locations[pi],
//                   index*n_samples_each_location+d_uniform_sampling_index[index],
//                   d_uniform_sampling[index*n_samples_each_location+d_uniform_sampling_index[index]],
//                   d_sum_weight[index]);
            while (d_uniform_sampling_index[index] < n_samples_each_location
            && d_uniform_sampling[index*n_samples_each_location+d_uniform_sampling_index[index]] < d_sum_weight[index]) {
//                printf("  while kernel location %d, n_distributions_each_location %d, pi %d, d_distribution_all_locations %f sum_weight %f "
//                       "d_uniform_sampling_index %d d_uniform_sampling[%d] = %f d_sample_index %d d_all_objects_index %d\n",
//                       index, n_distributions_each_location, pi, d_distribution_all_locations[pi], d_sum_weight[index],
//                       d_uniform_sampling_index[index],
//                       index*n_samples_each_location+d_uniform_sampling_index[index],
//                       d_uniform_sampling[index*n_samples_each_location+d_uniform_sampling_index[index]],
//                       d_sample_index[index*n_samples_each_location+d_uniform_sampling_index[index]],
//                       d_all_objects_index[pi]);
                d_sample_index[index*n_samples_each_location+d_uniform_sampling_index[index]] = d_all_objects_index[pi];
                d_uniform_sampling_index[index]++;
            }
            if (d_uniform_sampling_index[index] == n_samples_each_location) {
                return;
            }
        }
        __syncthreads();
    }
}

/*
 * This is GPU version of Random::roulette_sampling
 * d_n_samples size is n_locations
 * d_distributions size is n_locations*n_samples_each_location
 * all_objects size is n_locations*n_samples_each_location
 * d_sum_distribution size is n_locations
 * return size is n_locations*n_samples_each_location
 * */

template
TVector<Person*> GPU::Random::roulette_sampling<Person>(int, int,TVector<Person*>,ThrustTVectorDevice<double>,
                                                       ThrustTVectorDevice<double>&,bool);
template <class T>
TVector<T*> GPU::Random::roulette_sampling(int n_locations, int n_samples_each_location,
                                           TVector<T*> all_objects,
                                           ThrustTVectorDevice<double> d_distribution,
                                           ThrustTVectorDevice<double> &d_sum_distribution,
                                           bool is_shuffled){
    TVector<T*> samples(n_locations*n_samples_each_location, nullptr);
    double d_sum = thrust::reduce(thrust::device, d_sum_distribution.begin(), d_sum_distribution.end(), 0.0, thrust::plus<double>());
    if(d_sum == 0.0){
        return samples;
    }else if(d_sum == n_locations*(-1)){
        ThrustTVectorHost<double> h_sum(n_locations);
        for(int i = 0; i < n_locations; i++){
            int index_from = i*d_distribution.size();
            int index_to = index_from + d_distribution.size();
            h_sum[i] = thrust::reduce(thrust::device,
                                                d_distribution.begin() + index_from,
                                                d_distribution.begin() + index_to,
                                                0.0, thrust::plus<double>());
//            printf("GPU h_sum[%d] = %f\n",i,h_sum[i]);
        }
        d_sum_distribution = h_sum;
    }

    ThrustTVectorDevice<double> d_uniform_sampling(n_locations*n_samples_each_location);
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int n_blocks = (d_uniform_sampling.size() + n_threads + 1) / n_threads;
    random_uniform_kernel<<<n_blocks, n_threads>>>(d_states,
                                                   n_locations,
                                                   n_samples_each_location,
                                                   thrust::raw_pointer_cast(d_sum_distribution.data()),
                                                   thrust::raw_pointer_cast(d_uniform_sampling.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaPeekAtLastError());

//    thrust::copy(d_uniform_sampling.begin(), d_uniform_sampling.end(), std::ostream_iterator<double>(std::cout, "\n"));
//    printf("\n");

    for(int i = 0; i < n_locations; i++){
        int index_from = i*n_samples_each_location;
        int index_to = index_from + n_samples_each_location;
//        printf("location %d sort from %d to %d\n",i,index_from,index_to);
        thrust::sort(thrust::device,
                     d_uniform_sampling.begin() + index_from,
                     d_uniform_sampling.begin() + index_to);
    }
//    thrust::copy(d_uniform_sampling.begin(), d_uniform_sampling.end(), std::ostream_iterator<double>(std::cout, "\n"));
//    printf("\n");

    ThrustTVectorDevice<int> d_sample_index(n_locations*n_samples_each_location,0);
    ThrustTVectorDevice<double> d_sum_weight(n_locations,0.0);
    ThrustTVectorDevice<int> d_all_objects_index(n_locations*d_distribution.size(),0);
    thrust::sequence(thrust::device, d_all_objects_index.begin(), d_all_objects_index.end(), 0, 1);
    ThrustTVectorDevice<int> d_uniform_sampling_index(n_locations,0);
    roulette_sampling_kernel<<<n_blocks, n_threads>>>(n_locations,
                                                      n_samples_each_location,
                                                      d_distribution.size()/n_locations,
                                                      thrust::raw_pointer_cast(d_distribution.data()),
                                                      thrust::raw_pointer_cast(d_sum_weight.data()),
                                                      thrust::raw_pointer_cast(d_all_objects_index.data()),
                                                      thrust::raw_pointer_cast(d_uniform_sampling.data()),
                                                      thrust::raw_pointer_cast(d_uniform_sampling_index.data()),
                                                      thrust::raw_pointer_cast(d_sample_index.data()));
    cudaDeviceSynchronize();
    check_cuda_error(cudaPeekAtLastError());

//    thrust::copy(d_sample_index.begin(), d_sample_index.end(), std::ostream_iterator<int>(std::cout, "\n"));
//    printf("\n");

    if(is_shuffled){
        thrust::default_random_engine g;
        for(int i = 0; i < n_locations; i++){
            int index_from = i*n_samples_each_location;
            int index_to = index_from + n_samples_each_location;
//            printf("location %d shuffle from %d to %d\n",i,index_from,index_to);
            thrust::shuffle(thrust::device,
                            d_sample_index.begin() + index_from,
                            d_sample_index.begin() + index_to,
                            g);
        }
    }

//    thrust::copy(d_sample_index.begin(), d_sample_index.end(), std::ostream_iterator<double>(std::cout, "\n"));
//    printf("\n");

    for(int i = 0; i < d_sample_index.size(); i++){
        samples[i] = all_objects[d_sample_index[i]];
    }
    return samples;
}

template <class T>
TVector<T*> GPU::Random::multinomial_sampling(int n_locations, int n_samples_each_location,
                                           TVector<T*> all_objects,
                                           ThrustTVectorDevice<double> d_distribution,
                                           ThrustTVectorDevice<double> &d_sum_distribution,
                                           bool is_shuffled){

}


