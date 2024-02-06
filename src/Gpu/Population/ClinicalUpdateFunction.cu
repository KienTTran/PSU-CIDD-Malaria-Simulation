/* 
 * File:   ClinicalUpdateFunction.cpp
 * Author: Merlin
 * 
 * Created on July 29, 2013, 5:43 PM
 */

#include "ClinicalUpdateFunction.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

GPU::ClinicalUpdateFunction::ClinicalUpdateFunction(Model *model) : model_(model) {
}

GPU::ClinicalUpdateFunction::~ClinicalUpdateFunction() = default;

double GPU::ClinicalUpdateFunction::get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) {
  return model_->config()->parasite_density_level().log_parasite_density_asymptomatic;
}

