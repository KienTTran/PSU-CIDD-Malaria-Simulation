//
// Created by kient on 6/17/2023.
//

#ifndef MASS_RANDOM_CUH
#define MASS_RANDOM_CUH

#include <curand_kernel.h>
#include "Core/TypeDef.h"

namespace GPU{
    class Random;
}

class GPU::Random {
public:
    curandState *d_states;
    int n_threads;
    int n_blocks;
    ThrustTVectorDevice<curandState> d_location_curand_states;
public:
    Random();
    ~Random();

    void init(int n, long seed);
    void init_curand_states(ThrustTVectorDevice<curandState> &d_curand_states, int size, long seed);
    void free();
    void random_multinomial(int n_locations, int n_draws, ThrustTVectorDevice<int> N, ThrustTVectorDevice<double> &d_p, ThrustTVectorDevice<unsigned int> &d_n);
    void sampling_multinomial(int n_locations, int K, ThrustTVectorDevice<int> N, ThrustTVectorDevice<double> &h_p, ThrustTVectorHost<unsigned int> &h_samples);

    int randomPoisson(int n, float *A, int *B, int nthreads, unsigned long long seed, unsigned long long offset);
    int randomBinomial(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out, unsigned long long seed, unsigned long long offset);
    int randomGamma(int nrows, int ncols, float *A, int atype, float *B, int btype, float *Out, unsigned long long seed, unsigned long long offset);
};


#endif //MASS_RANDOM_CUH
