/* 
 * File:   ImmunityClearanceUpdateFunction.cpp
 * Author: Merlin
 * 
 * Created on July 29, 2013, 5:49 PM
 */

#include "ImmunityClearanceUpdateFunction.cuh"
#include "ClonalParasitePopulation.cuh"
#include "SingleHostClonalParasitePopulations.cuh"
#include "Person.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"

__device__ __host__ GPU::ImmunityClearanceUpdateFunction::ImmunityClearanceUpdateFunction(Model *model) : model_(model) {}

GPU::ImmunityClearanceUpdateFunction::~ImmunityClearanceUpdateFunction() = default;

__host__ double GPU::ImmunityClearanceUpdateFunction::get_current_parasite_density(GPU::ClonalParasitePopulation *parasite,
                                                                                              int duration) {

  auto *p = parasite->parasite_population()->person();
  return p->immune_system()->get_parasite_size_after_t_days(duration, parasite->last_update_log10_parasite_density(),
                                                            parasite->genotype()->daily_fitness_multiple_infection);
}

__device__ double GPU::ImmunityClearanceUpdateFunction::get_current_parasite_density_gpu(ParasiteDensityLevel h_parasite_density_level,
                                                                                         double person_immune_parasite_size_after_t_days,
                                                                                         int duration){
    printf("GPU::ImmunityClearanceUpdateFunction::get_current_parasite_density_gpu %f\n",person_immune_parasite_size_after_t_days);
    return person_immune_parasite_size_after_t_days;
//  return p->immune_system()->get_parasite_size_after_t_days(duration, parasite->last_update_log10_parasite_density(),
//                                                            parasite->genotype()->daily_fitness_multiple_infection);
}
