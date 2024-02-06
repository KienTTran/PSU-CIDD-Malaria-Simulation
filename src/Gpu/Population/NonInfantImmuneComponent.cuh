/* 
 * File:   NormalImmuneComponent.h
 * Author: nguyentran
 *
 * Created on May 28, 2013, 11:06 AM
 */

#ifndef NONINFANTIMMUNECOMPONENT_CUH
#define    NONINFANTIMMUNECOMPONENT_CUH

#include "ImmuneComponent.cuh"
#include "Core/TypeDef.h"

namespace GPU{
    class NonInfantImmuneComponent;
    class ImmuneSystem;
}

class GPU::NonInfantImmuneComponent : public GPU::ImmuneComponent {
 public:
    NonInfantImmuneComponent(GPU::ImmuneSystem *immune_system = nullptr);

  // NonInfantImmuneComponent(const NonInfantImmuneComponent& orig);
  virtual ~NonInfantImmuneComponent();

    virtual double get_decay_rate(ImmuneSystemInformation h_immune_system_information,const int &age = 0) const;

    virtual double get_acquire_rate(ImmuneSystemInformation h_immune_system_information,int latest_update_time,const int &age = 0) const;

    virtual int type() const {
        return 2;
    }
public:
    ImmuneSystemInformation h_immune_system_information;
 private:

};

#endif    /* NONINFANTIMMUNECOMPONENT_CUH */

