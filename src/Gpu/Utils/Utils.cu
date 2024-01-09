//
// Created by kient on 6/17/2023.
//

#include "Utils.cuh"
#include "Population/Person.h"

#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/sequence.h>
#include <thrust/remove.h>
#include "Model.h"
#include "Core/Config/Config.h"

GPU::Utils::Utils() {
}

GPU::Utils::~Utils() {
}

void GPU::Utils::init(){
}

struct CheckKeyT2
{
    template<typename T,typename T2>
    __host__ __device__
    bool operator()(thrust::tuple<T,T2> t0, thrust::tuple<T,T2> t1)
    {
        if(thrust::get<0>(t0) == thrust::get<0>(t1)) return true;
        else return false;
    }
};

struct SumValueT2
{
    template<typename T,typename T2>
    __host__ __device__
    thrust::tuple<T, T2> operator()(thrust::tuple<T,T2> t0, thrust::tuple<T,T2> t1)
    {
        return thrust::make_tuple(thrust::get<0>(t1),thrust::get<1>(t0) + thrust::get<1>(t1));
    }
};

struct IsValueZeroT2
{
    template<typename T,typename T2>
    __host__ __device__
    bool operator()(thrust::tuple<T,T2> t)
    {
        return thrust::get<1>(t) == 0;
    }
};

/*
 * This function is used to sum the values of an array grouped by value in another array
 * Two arrays must have the same size and size parameter is the length of output array
 * Output is a tuple with 2 elements: key and sum value of that key
 * For example:
 * input_keys = [0, 1, 2, 1, 3, 0] - location
 * input_values = [1, 2, 3, 4, 5, 6] - age
 * output = [(0,7), (1,6), (2,3), (3,5)]
*/
template //https://stackoverflow.com/a/51606460/9187675 - to call template functions from other files
ThrustTuple2Vector<int,int> GPU::Utils::sum_value_by_1key<int,int>(TVector<int>, TVector<int>, int);
template
ThrustTuple2Vector<double,double> GPU::Utils::sum_value_by_1key<double,double>(TVector<double>, TVector<double>, int);
template
ThrustTuple2Vector<int,double> GPU::Utils::sum_value_by_1key<int,double>(TVector<int>, TVector<double>, int);
template
ThrustTuple2Vector<double,int> GPU::Utils::sum_value_by_1key<double,int>(TVector<double>, TVector<int>, int);
template<typename T,typename T2>
ThrustTuple2Vector<T,T2> GPU::Utils::sum_value_by_1key(TVector<T> input_keys, TVector<T2> input_values, int size){

    ThrustTVectorDevice<T> device_keys = input_keys;
    ThrustTVectorDevice<T2> device_values = input_values;

    thrust::sort_by_key(thrust::device, device_keys.begin(), device_keys.end(), device_values.begin(), thrust::less<T>());

    auto begin = thrust::make_zip_iterator(thrust::make_tuple(device_keys.begin(), device_values.begin()));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(device_keys.end(), device_values.end()));

    ThrustTuple2VectorDevice<T,T2> device_output_values(device_values.size());

    auto result = thrust::reduce_by_key(thrust::device,
                                        begin,
                                        end,
                                        begin,
                                        thrust::make_discard_iterator(),
                                        device_output_values.begin(),
                                        CheckKeyT2(),
                                        SumValueT2());
    int output_length = result.second - device_output_values.begin();
    ThrustTuple2Vector<T,T2> host_output_values(output_length);
    thrust::copy(device_output_values.begin(), device_output_values.begin() + output_length, host_output_values.begin());
    thrust::remove_if(host_output_values.begin(), host_output_values.end(), IsValueZeroT2());
    return host_output_values;
}

template<typename T>
__global__ void fill_missing_indices(thrust::tuple<T,int>* device_output_values, T* output, int size){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < size; i += stride){
        output[thrust::get<0>(device_output_values[i])] = thrust::get<1>(device_output_values[i]);
    }
}

/*
 * This function is used to count occurrence of each key in an array
 * Input is key array with size parameter is the length of output array
 * For example:
 * input_keys = [0, 1, 2, 1, 3, 0]
 * size = 4
 * output = [2, 2, 1, 1]
 * Index which is not in the key array will be count as 0
*/
template
TVector<int> GPU::Utils::count_by_1key<int>(TVector<int>, int);
template<typename T>
TVector<T> GPU::Utils::count_by_1key(TVector<T> input_keys, int size){
    ThrustTVectorDevice<T> device_keys = input_keys;
    thrust::sort(thrust::device, device_keys.begin(), device_keys.end(), thrust::less<T>());
    auto begin = thrust::make_zip_iterator(thrust::make_tuple(device_keys.begin(), device_keys.begin()));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(device_keys.end(), device_keys.end()));

    ThrustTuple2VectorDevice<T,int> device_output_temp(input_keys.size());

    auto result = thrust::reduce_by_key(thrust::device,
                                        begin,
                                        end,//https://stackoverflow.com/questions/34250322/count-reduction-using-thrust
                                        thrust::make_zip_iterator(thrust::make_tuple(device_keys.begin(),thrust::make_constant_iterator(1))),
                                        thrust::make_discard_iterator(),
                                        device_output_temp.begin(),
                                        CheckKeyT2(),
                                        SumValueT2());
    ThrustTVectorDevice<T> device_output_values(size,0);
    TVector<T> host_output_values(size,0);
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int n_blocks = (size + n_threads - 1) / n_threads;
    fill_missing_indices<T><<<n_blocks,n_threads>>>(thrust::raw_pointer_cast(device_output_temp.data()),
                                                      thrust::raw_pointer_cast(device_output_values.data()),
                                                      size);
    cudaDeviceSynchronize();
    check_cuda_error(cudaGetLastError());
    thrust::copy(device_output_values.begin(), device_output_values.end(), host_output_values.begin());
    return host_output_values;
}

/*
 * Get 3 vectors of from_location, to_location and moving_level for removing zero from
 * d_circulation_indices and d_n_circulations_all_loc_ml
 * */
__global__ void extract_locations_and_moving_levels_kernel(int n_locations,
                                                    int n_moving_levels,
                                                    thrust::tuple<int,int> *d_circulation_indices,
                                                    unsigned int *d_n_circulations_all_loc_ml,
                                                    int* d_from_indices,
                                                    int* d_to_indices,
                                                    int* d_moving_levels){
    int thread_index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations*n_moving_levels; index += stride) {
        if(d_n_circulations_all_loc_ml[index] == 0) return;
        int circulate_index = index / n_moving_levels;
        d_from_indices[index] = thrust::get<0>(d_circulation_indices[circulate_index]);
        d_to_indices[index] = thrust::get<1>(d_circulation_indices[circulate_index]);
        d_moving_levels[index] = index % n_moving_levels;
    }
}
