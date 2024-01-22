/*
 * File:   BloodParasite.cpp
 * Author: Merlin
 *
 * Created on July 11, 2013, 2:21 PM
 */

#include "ClonalParasitePopulation.cuh"

#include <cmath>

#include "Core/Config/Config.h"
#include "Helpers/NumberHelpers.h"
#include "Model.h"
#include "Therapies/Therapy.h"
#include "Gpu/Population/Person.cuh"
#include "ClinicalUpdateFunction.cuh"
#include "ImmuneSystem.cuh"
#include "Gpu/Utils/Utils.cuh"

const double GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY = -1000.0;

GPU::ClonalParasitePopulation::ClonalParasitePopulation(GPU::Genotype *genotype)
    : last_update_log10_parasite_density_(-1000.0),
      gametocyte_level_(0.0),
      first_date_in_blood_(-1),
      parasite_population_(nullptr),
      genotype_(genotype),
      update_function_(nullptr) {
}

GPU::ClonalParasitePopulation::~ClonalParasitePopulation() = default;

__host__ void GPU::ClonalParasitePopulation::set_gpu_update_function(GPU::ParasiteDensityUpdateFunction *h_function) {
}

__host__ double GPU::ClonalParasitePopulation::get_current_parasite_density(const int &current_time) {
//  printf("GPU::ClonalParasitePopulation::get_current_parasite_density %f\n",last_update_log10_parasite_density_);
  const auto duration = current_time - parasite_population_->person()->latest_update_time();
  if (duration == 0) {
      printf("duration == 0\n");
    return last_update_log10_parasite_density_;
  }

  if (update_function_ == nullptr) {
      printf("update_function_ == nullptr\n");
    //        std::cout << "hello" << std::endl;
    return last_update_log10_parasite_density_;
  }

  return update_function_->get_current_parasite_density(this, duration);
}

__device__ __host__ double GPU::ClonalParasitePopulation::get_log10_infectious_density() const {
  if (NumberHelpers::is_equal(last_update_log10_parasite_density_, -1000.0)
      || NumberHelpers::is_equal(gametocyte_level_, 0.0))
    return -1000.0;

  return last_update_log10_parasite_density_ + log10(gametocyte_level_);
}

double GPU::ClonalParasitePopulation::last_update_log10_parasite_density() const {
  return last_update_log10_parasite_density_;
}

void GPU::ClonalParasitePopulation::set_last_update_log10_parasite_density(const double &value) {
  last_update_log10_parasite_density_ = value;
}

double GPU::ClonalParasitePopulation::gametocyte_level() const {
  return gametocyte_level_;
}

void GPU::ClonalParasitePopulation::set_gametocyte_level(const double &value) {
  if (NumberHelpers::is_enot_qual(gametocyte_level_, value)) {
    gametocyte_level_ = value;
  }
}

GPU::Genotype *GPU::ClonalParasitePopulation::genotype() const {
  return genotype_;
}

void GPU::ClonalParasitePopulation::set_genotype(GPU::Genotype *value) {
  if (genotype_ != value) {
    genotype_ = value;
  }
}

bool GPU::ClonalParasitePopulation::resist_to(const int &drug_id) const {
  return genotype_->resist_to(Model::CONFIG->drug_db()->at(drug_id));
}

void GPU::ClonalParasitePopulation::update() {
  set_last_update_log10_parasite_density(get_current_parasite_density(Model::GPU_SCHEDULER->current_time()));
}

void GPU::ClonalParasitePopulation::perform_drug_action(const double &percent_parasite_remove) {
  double newSize = last_update_log10_parasite_density_;
  if (percent_parasite_remove > 1) {
    newSize = Model::CONFIG->parasite_density_level().log_parasite_density_cured;
  } else {
    newSize += log10(1 - percent_parasite_remove);
  }

  if (newSize < Model::CONFIG->parasite_density_level().log_parasite_density_cured) {
    newSize = Model::CONFIG->parasite_density_level().log_parasite_density_cured;
  }

  //    std::cout << Model::GPU_SCHEDULER->current_time() << "\t" <<parasite_population()->person() << "\t"  <<
  //    percent_parasite_remove << "\t"<<last_update_log10_parasite_density_ << "\t" <<newSize << std::endl;
  set_last_update_log10_parasite_density(newSize);
}

__device__ double GPU::ClonalParasitePopulation::get_current_parasite_density_gpu(GPU::ParasiteDensityUpdateFunction* d_update_function,
                                                                                  GPU::Genotype* d_genotype,
                                                                                  GPU::ImmuneSystem* d_immune_system,
                                                                                  ParasiteDensityLevel h_parasite_density_level,
                                                                                  ImmuneSystemInformation* d_immune_system_information,
                                                                                  const int current_time,
                                                                                  const int latest_update_time) {
    const auto duration = current_time - latest_update_time;
    if (duration == 0) {
        printf("duration == 0\n");
        return last_update_log10_parasite_density_;
    }

//    printf("GPU::ClonalParasitePopulation::get_current_parasite_density_gpu %d %d %d %f\n",
//           current_time,latest_update_time,duration,last_update_log10_parasite_density_);

    return 5.6;
//    if (d_update_function == nullptr) {
//        printf("d_update_function ==  nullptr\n");
//        return last_update_log10_parasite_density_;
//    }
//
//    double parsite_size = d_immune_system->get_parasite_size_after_t_days_gpu(d_immune_system_information,duration,
//                                                                              last_update_log10_parasite_density_,
//                                                                              d_genotype->daily_fitness_multiple_infection);
//    return d_update_function->get_current_parasite_density_gpu(h_parasite_density_level,
//                                                               parsite_size,
//                                                               current_time);
}

__device__ __host__ double GPU::ClonalParasitePopulation::test2(){
    test_ = 99.9;
    return 45.45;
}
