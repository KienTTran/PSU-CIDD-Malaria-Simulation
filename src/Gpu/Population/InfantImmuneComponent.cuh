/* 
 * File:   InfantImmuneComponent.h
 * Author: nguyentran
 *
 * Created on May 28, 2013, 11:06 AM
 */

#ifndef INFANTIMMUNECOMPONENT_CUH
#define    INFANTIMMUNECOMPONENT_CUH

#include "ImmuneComponent.cuh"

namespace GPU {
  class InfantImmuneComponent;
}

class GPU::InfantImmuneComponent : public GPU::ImmuneComponent {
 public:
  explicit __device__ __host__ InfantImmuneComponent(GPU::ImmuneSystem *immune_system = nullptr);

  // InfantImmuneComponent(const InfantImmuneComponent& orig);
  virtual ~InfantImmuneComponent();

    __device__ __host__ double get_decay_rate(ImmuneSystemInformation h_immune_system_information,const int &age) const override;

    __device__ __host__ double get_acquire_rate(ImmuneSystemInformation h_immune_system_information,int latest_update_time,const int &age) const override;

    __device__ __host__ double get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) override;

};

#endif    /* INFANTIMMUNECOMPONENT_CUH */
