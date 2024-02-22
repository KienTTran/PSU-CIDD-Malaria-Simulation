/* 
 * File:   NormalImmuneComponent.cu
 * Author: nguyentran
 * 
 * Created on May 28, 2013, 11:06 AM
 */

#include "NonInfantImmuneComponent.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

GPU::NonInfantImmuneComponent::NonInfantImmuneComponent(GPU::ImmuneSystem *immune_system) : GPU::ImmuneComponent(immune_system) {
}

GPU::NonInfantImmuneComponent::~NonInfantImmuneComponent() {
}

double GPU::NonInfantImmuneComponent::get_acquire_rate(ImmuneSystemInformation h_immune_system_information,int latest_update_time,const int &age) const {
  //    return FastImmuneComponent::acquireRate;

  return (age > 80) ? h_immune_system_information.acquire_rate_by_age[80]
                    : h_immune_system_information.acquire_rate_by_age[age];

}

double GPU::NonInfantImmuneComponent::get_decay_rate(ImmuneSystemInformation h_immune_system_information,const int &age) const {
  return h_immune_system_information.decay_rate;
}