//
// Created by kient on 6/17/2023.
//

#ifndef RANDOM_CUH
#define RANDOM_CUH

#include <curand_kernel.h>
#include "Core/TypeDef.h"

namespace GPU{
    class Random;
}

class GPU::Random {
public:
    curandState *d_states;
public:
    Random();
    ~Random();

    void init(int n, unsigned long seed, int n_threads = -1,bool debug = false);
    void init_curand_states(ThrustTVectorDevice<curandState> &d_curand_states, int size, long seed);
    void free();
    void random_multinomial(int n_locations, int n_samples_each_location,
                                         ThrustTVectorDevice<int> d_n_trials,
                                         ThrustTVectorDevice<double> d_distributions,
                                         ThrustTVectorDevice<unsigned int> &d_samples);
    template <class T>
    TVector<T*> roulette_sampling(int n_locations, int n_samples_each_location,
                                  ThrustTVectorDevice<double> d_distribution_all_locations,
                                  TVector<T*> all_objects,
                                  ThrustTVectorDevice<double> &d_sum_distribution_all_locations,
                                  bool is_shuffled);
    template <class T>
    TVector<T*> multinomial_sampling(int n_locations, int n_samples_each_location,
                                     ThrustTVectorDevice<double> d_distribution_all_locations,
                                     TVector<T*> all_objects,
                                     ThrustTVectorDevice<double> &d_sum_distribution_all_locations,
                                     bool is_shuffled);
    ThrustTVectorDevice<double> random_uniform_double_min_max(int size, double from, double to);
    ThrustTVectorDevice<int> random_uniform_int_min_max(int size, int from, int to);
};


#endif //RANDOM_CUH
