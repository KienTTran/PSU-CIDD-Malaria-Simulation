/* 
 * File:   ImmuneSystem.cpp
 * Author: nguyentran
 * 
 * Created on May 27, 2013, 11:56 AM
 */

#include "ImmuneSystem.cuh"
#include "ImmuneComponent.cuh"
#include "Person.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include <cmath>
#include "Helpers/ObjectHelpers.h"

__device__ __host__ GPU::ImmuneSystem::ImmuneSystem(GPU::Person *p) : person_(p), increase_(false), immune_component_(nullptr) {
  //    immune_components_ = new ImmuneComponentPtrVector();

}

GPU::ImmuneSystem::~ImmuneSystem() {

  if (immune_component_!=nullptr) {
    ObjectHelpers::delete_pointer<ImmuneComponent>(immune_component_);
  }
  assert(immune_component_==nullptr);
  person_ = nullptr;
}

GPU::ImmuneComponent *GPU::ImmuneSystem::immune_component() const {
  return immune_component_;
}

void GPU::ImmuneSystem::set_immune_component(ImmuneComponent *value) {
  if (immune_component_!=value) {
    if (immune_component_!=nullptr) {
      ObjectHelpers::delete_pointer<ImmuneComponent>(immune_component_);
    }

    immune_component_ = value;
    immune_component_->set_immune_system(this);
  }
}

void GPU::ImmuneSystem::draw_random_immune(double value) {
  immune_component_->draw_random_immune(value);
}

__device__ __host__ double GPU::ImmuneSystem::get_lastest_immune_value() const {
  return immune_component_->latest_value();
}

void GPU::ImmuneSystem::set_latest_immune_value(double value) {
  immune_component_->set_latest_value(value);
}

__device__ __host__ double GPU::ImmuneSystem::get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) const {
  return immune_component_->get_current_value(h_immune_system_information,latest_update_time,current_time);
}

__host__ double GPU::ImmuneSystem::get_parasite_size_after_t_days(const int &duration, const double &originalSize,
                                                    const double &fitness) const {

  const auto last_immune_level = get_lastest_immune_value();
  const auto temp = Model::CONFIG->immune_system_information().c_max*(1 - last_immune_level) + Model::CONFIG
      ->immune_system_information()
      .c_min*
      last_immune_level;

  const auto value = originalSize + duration*(log10(temp) + log10(fitness));
  return value;

}

__device__ __host__ double GPU::ImmuneSystem::get_parasite_size_after_t_days_gpu(ImmuneSystemInformation* immunity_system_info,
                                                                        const int &duration, const double &originalSize,
                                                                        const double &fitness) const {
    printf("GPU::ImmuneSystem::get_parasite_size_after_t_days_gpu %d %f %f %f %f %f\n",
           duration,immunity_system_info->c_max,immunity_system_info->c_min,originalSize,fitness,get_lastest_immune_value());
    const auto last_immune_level = get_lastest_immune_value();
    const auto temp = immunity_system_info->c_max*(1 - last_immune_level) + immunity_system_info->c_min*last_immune_level;
    const auto value = originalSize + duration*(log10(temp) + log10(fitness));
    return value;
}


__device__ __host__ double GPU::ImmuneSystem::get_clinical_progression_probability(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) const {
  const auto immune = get_current_value(h_immune_system_information,latest_update_time,current_time);
  //    double PClinical = (isf.min_clinical_probability - isf.max_clinical_probability) * pow(immune, isf.immune_effect_on_progression_to_clinical) + isf.max_clinical_probability;

  //    const double p_m = 0.99;

  const double mid_point = 0.4;
  const auto p_clinical = h_immune_system_information.max_clinical_probability/(1 + pow((immune/mid_point),
                                                                h_immune_system_information.immune_effect_on_progression_to_clinical));

  //    std::cout << immune << "\t" << PClinical<< std::endl;
  return p_clinical;
}

__device__ __host__ void GPU::ImmuneSystem::update(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) {
  immune_component_->update(h_immune_system_information,latest_update_time,current_time);
}
