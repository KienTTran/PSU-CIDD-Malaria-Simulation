/* 
 * File:   ImmuneComponent.cpp
 * Author: nguyentran
 * 
 * Created on May 27, 2013, 12:44 PM
 */

#include "ImmuneComponent.cuh"
#include "ImmuneSystem.cuh"
#include "Person.cuh"
#include "Gpu/Population/Population.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include <cmath>

__device__ __host__ GPU::ImmuneComponent::ImmuneComponent(GPU::ImmuneSystem *immune_system) : immune_system_(immune_system), latest_value_(0.0) {}

GPU::ImmuneComponent::~ImmuneComponent() {
  immune_system_ = nullptr;
}

__device__ __host__ double GPU::ImmuneComponent::get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) {
  auto temp = 0.0;
  if (immune_system_!=nullptr) {
    if (immune_system_->person()!=nullptr) {
      const auto duration = current_time - latest_update_time;

      const auto age = immune_system_->person()->age();
      if (immune_system_->increase()) {
        //increase I(t) = 1 - (1-I0)e^(-b1*t)

        temp = 1 - (1 - latest_value_)*exp(-get_acquire_rate(h_immune_system_information,age)*duration);

        //        temp = lastImmuneLevel;
        //        double b1 = GetAcquireRate(immuneSystem->person->age);
        //        for (int i = 0; i < duration; i++) {
        //            temp+= b1*(1-temp);
        //        }

      } else {
        //decrease I(t) = I0 * e ^ (-b2*t);
        temp = latest_value_*exp(-get_decay_rate(h_immune_system_information,age)*duration);
        temp = (temp < 0.00001) ? 0.0 : temp;
      }

    }
  }
  return temp;
}

__device__ __host__ void GPU::ImmuneComponent::update(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) {
  latest_value_ = get_current_value(h_immune_system_information,latest_update_time,current_time);

}

/*
 * Assign value instead of draw random here
 * */
__device__ __host__ void GPU::ImmuneComponent::draw_random_immune(double value) {
  //TODO: Double check this is used or not
  latest_value_ = value;//Model::RANDOM->random_beta(h_immune_system_information.alpha_immune, h_immune_system_information.beta_immune);
}
