/* 
 * File:   ImmuneSystem.h
 * Author: nguyentran
 *
 * Created on May 27, 2013, 11:56 AM
 */

#ifndef IMMUNESYSTEM_CUH
#define    IMMUNESYSTEM_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

class Model;

namespace GPU{
    class Person;
    class ImmuneSystem;
    class ImmuneComponent;
}

class Config;

//typedef std::vector<ImmuneComponent*> ImmuneComponentPtrVector;

class GPU::ImmuneSystem {

 POINTER_PROPERTY(GPU::Person, person)

 VIRTUAL_PROPERTY_REF(bool, increase)
  //    POINTER_PROPERTY_REF(ImmuneComponentPtrVector, immune_components);
 POINTER_PROPERTY_HEADER(GPU::ImmuneComponent, immune_component)

 public:
  explicit __device__ __host__ ImmuneSystem(GPU::Person *p = nullptr);

  virtual ~ImmuneSystem();

  //    virtual void clear();

  virtual void draw_random_immune(double value);

  virtual __device__ __host__ void update(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time);

  virtual __device__ __host__ double get_lastest_immune_value() const;

  virtual void set_latest_immune_value(double value);

  virtual __device__ __host__ double get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) const;

  virtual __host__ double
  get_parasite_size_after_t_days(const int &duration, const double &originalSize, const double &fitness) const;

  virtual __device__ __host__ double
  get_parasite_size_after_t_days_gpu(ImmuneSystemInformation* h_immune_system_information,
                                     const int &duration, const double &originalSize,
                                     const double &fitness) const;

  virtual __device__ __host__ double get_clinical_progression_probability(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time) const;

public:
    double test_ = 3.0;
    virtual __device__ __host__ double test(ImmuneSystemInformation* d_immune_system_information){
      test_ = 7.0;
      return test_;
    }
};

#endif    /* IMMUNESYSTEM_CUH */
