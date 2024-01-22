/* 
 * File:   ClinicalUpdateFunction.cpp
 * Author: Merlin
 * 
 * Created on July 29, 2013, 5:43 PM
 */

#include "ClinicalUpdateFunction.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

__device__ __host__ GPU::ClinicalUpdateFunction::ClinicalUpdateFunction(Model *model) : model_(model) {
}

GPU::ClinicalUpdateFunction::~ClinicalUpdateFunction() = default;

__host__ double GPU::ClinicalUpdateFunction::get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) {
  return model_->config()->parasite_density_level().log_parasite_density_asymptomatic;
}

__device__ double GPU::ClinicalUpdateFunction::get_current_parasite_density_gpu(ParasiteDensityLevel h_parasite_density_level,
                                                                                double person_immune_parasite_size_after_t_days,
                                                                                int duration) {
//  printf("GPU::ClinicalUpdateFunction::get_current_parasite_density_gpu %f\n",h_parasite_density_level.log_parasite_density_asymptomatic);
  return h_parasite_density_level.log_parasite_density_asymptomatic;
}
