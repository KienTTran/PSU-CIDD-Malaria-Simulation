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
#include "Gpu/Therapies/DrugDatabase.cuh"
#include "Gpu/Therapies/Therapy.cuh"
#include "Gpu/Population/Person.cuh"
#include "ClinicalUpdateFunction.cuh"
#include "ImmuneSystem.cuh"
#include "Gpu/Utils/Utils.cuh"
#include "Population.cuh"
#include "Properties/PersonIndexGPU.cuh"

GPU::ClonalParasitePopulation::ClonalParasitePopulation(GPU::Genotype *genotype)
    : last_update_log10_parasite_density_(LOG_ZERO_PARASITE_DENSITY),
      gametocyte_level_(0.0),
      first_date_in_blood_(-1),
      parasite_population_(nullptr),
      genotype_(genotype),
      id_(UniqueId::get_instance().get_uid()),
      index_(-1),
      person_(nullptr){
}

GPU::ClonalParasitePopulation::~ClonalParasitePopulation() = default;


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

  if(person_->index() >= 1040 && person_->index() < 1010){
    printf("%d GPU::ClonalParasitePopulation::get_current_parasite_density %d %d %d %d %d %d %f\n",
           person_->index(),
           index_,
           parasite_population_->person()->latest_update_time(),
           person_->person_index_gpu->h_person_update_info()[person_->index()].person_latest_update_time,
           current_time,
           duration,
           update_function_->type(),
           update_function_->get_current_parasite_density(this, duration));
  }
  return update_function_->get_current_parasite_density(this, duration);
}

__device__ __host__ double GPU::ClonalParasitePopulation::get_log10_infectious_density() const {
  if (NumberHelpers::is_equal(last_update_log10_parasite_density_, LOG_ZERO_PARASITE_DENSITY)
      || NumberHelpers::is_equal(gametocyte_level_, 0.0))
    return LOG_ZERO_PARASITE_DENSITY;

  return last_update_log10_parasite_density_ + log10(gametocyte_level_);
}

double GPU::ClonalParasitePopulation::last_update_log10_parasite_density() const {
  return last_update_log10_parasite_density_;
}

void GPU::ClonalParasitePopulation::set_last_update_log10_parasite_density(const double &value) {
  last_update_log10_parasite_density_ = value;
  person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_last_update_log10_parasite_density[index_]
  = last_update_log10_parasite_density_;
}

double GPU::ClonalParasitePopulation::gametocyte_level() const {
  return gametocyte_level_;
}

void GPU::ClonalParasitePopulation::set_gametocyte_level(const double &value) {
  if (NumberHelpers::is_enot_qual(gametocyte_level_, value)) {
    gametocyte_level_ = value;
    person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_gametocyte_level[index_] = gametocyte_level_;
  }
}

GPU::Genotype *GPU::ClonalParasitePopulation::genotype() const {
  return genotype_;
}

void GPU::ClonalParasitePopulation::set_genotype(GPU::Genotype *value) {
  if (genotype_ != value) {
    genotype_ = value;
    /* Convert string to char[] with NULL terminated */
    std::copy( genotype_->aa_sequence.begin(),
               genotype_->aa_sequence.end(),
               person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_genotype[index_]);
    person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_genotype[index_][MAX_GENOTYPE_LOCUS] = '\0';
    person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_genotype_fitness_multiple_infection[index_]
    = genotype_->daily_fitness_multiple_infection;
  }
}

bool GPU::ClonalParasitePopulation::resist_to(const int &drug_id) const {
  return genotype_->resist_to(Model::CONFIG->gpu_drug_db()->at(drug_id));
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

void GPU::ClonalParasitePopulation::set_update_function(GPU::ParasiteDensityUpdateFunction *value) {
  update_function_ = value;
  person_->person_index_gpu->h_person_update_info()[person_->index()].parasite_update_function_type[index_] = value->type();
}

GPU::ParasiteDensityUpdateFunction *GPU::ClonalParasitePopulation::update_function() const {
  return update_function_;
}
