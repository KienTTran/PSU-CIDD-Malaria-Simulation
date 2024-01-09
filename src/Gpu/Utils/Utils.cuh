//
// Created by kient on 6/17/2023.
//

#ifndef MASS_UTILS_CUH
#define MASS_UTILS_CUH

#include <iostream>
#include <curand_kernel.h>
#include "Core/TypeDef.h"
#include "Helpers/UniqueId.hxx"

#define check_cuda_error(ans) { gpu_assert((ans), __FILE__, __LINE__); }
inline void gpu_assert(cudaError_t code, const char* file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        printf("check_cuda_error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
namespace GPU{
    class Utils;
}

class GPU::Utils {
public:
    Utils();
    ~Utils();
    void init();
public:
    template<typename T,typename T2>
    ThrustTuple2Vector<T,T2> sum_value_by_1key(TVector<T> input_keys, TVector<T2> input_values, int size);
    template<typename T>
    TVector<T> count_by_1key(TVector<T> input_keys, int size);
public:
};
//
////From https://github.com/NVIDIA/thrust/blob/master/examples/strided_range.cu
//#include <thrust/iterator/counting_iterator.h>
//#include <thrust/iterator/transform_iterator.h>
//#include <thrust/iterator/permutation_iterator.h>
//#include <thrust/functional.h>
//#include <thrust/fill.h>
//#include <thrust/device_vector.h>
//#include <thrust/copy.h>
//#include <iostream>
//
//// this example illustrates how to make strided access to a range of values
//// examples:
////   strided_range([0, 1, 2, 3, 4, 5, 6], 1) -> [0, 1, 2, 3, 4, 5, 6]
////   strided_range([0, 1, 2, 3, 4, 5, 6], 2) -> [0, 2, 4, 6]
////   strided_range([0, 1, 2, 3, 4, 5, 6], 3) -> [0, 3, 6]
////   ...
//
//template <typename Iterator>
//class strided_range
//{
//public:
//
//    typedef typename thrust::iterator_difference<Iterator>::type difference_type;
//
//    struct stride_functor : public thrust::unary_function<difference_type,difference_type>
//    {
//        difference_type stride;
//
//        stride_functor(difference_type stride)
//                : stride(stride) {}
//
//        __host__ __device__
//        difference_type operator()(const difference_type& i) const
//        {
//            return stride * i;
//        }
//    };
//
//    typedef typename thrust::counting_iterator<difference_type>                   CountingIterator;
//    typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
//    typedef typename thrust::permutation_iterator<Iterator,TransformIterator>     PermutationIterator;
//
//    // type of the strided_range iterator
//    typedef PermutationIterator iterator;
//
//    // construct strided_range for the range [first,last)
//    strided_range(Iterator first, Iterator last, difference_type stride)
//            : first(first), last(last), stride(stride) {}
//
//    iterator begin(void) const
//    {
//        return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
//    }
//
//    iterator end(void) const
//    {
//        return begin() + ((last - first) + (stride - 1)) / stride;
//    }
//
//protected:
//    Iterator first;
//    Iterator last;
//    difference_type stride;
//};


#endif //MASS_UTILS_CUH
