/* 
 * File:   InfantImmuneComponent.cpp
 * Author: nguyentran
 * 
 * Created on May 28, 2013, 11:06 AM
 */

#include "InfantImmuneComponent.cuh"
#include "Model.h"
#include "Gpu/Population/Person.cuh"
#include "ImmuneSystem.cuh"
#include <cmath>


__device__ __host__ GPU::InfantImmuneComponent::InfantImmuneComponent(GPU::ImmuneSystem *immune_system) : GPU::ImmuneComponent(immune_system) {}

GPU::InfantImmuneComponent::~InfantImmuneComponent() = default;

__device__ __host__ double GPU::InfantImmuneComponent::get_acquire_rate(ImmuneSystemInformation h_immune_system_information,int latest_update_time,const int &age) const {
  return 0;
}

__device__ __host__ double GPU::InfantImmuneComponent::get_decay_rate(ImmuneSystemInformation h_immune_system_information,const int &age) const {
  return 0.0315;
}

__device__ __host__ double GPU::InfantImmuneComponent::get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) {
  auto temp = 0.0;
  if (latest_update_time > 0) {
    const auto duration = current_time - latest_update_time;
    //decrease I(t) = I0 * e ^ (-b2*t);
    temp = latest_value()*exp(-get_decay_rate(h_immune_system_information,0)*duration);
  }
  return temp;
}
