/* 
 * File:   ImmuneComponent.h
 * Author: nguyentran
 *
 * Created on May 27, 2013, 12:44 PM
 */

#ifndef IMMUNECOMPONENT_CUH
#define    IMMUNECOMPONENT_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

//#include "ObjectPool.h"
namespace GPU{
    class ImmuneComponent;
    class ImmuneSystem;
}

class Model;

class GPU::ImmuneComponent {
 POINTER_PROPERTY(GPU::ImmuneSystem, immune_system)

 PROPERTY_REF(double, latest_value)

 public:
  explicit ImmuneComponent(GPU::ImmuneSystem *immune_system = nullptr);

  //    ImmuneComponent(const ImmuneComponent& orig);
  virtual ~ImmuneComponent();

  void update(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time);

  virtual void draw_random_immune(double value);

  virtual double get_current_value(ImmuneSystemInformation h_immune_system_information,int latest_update_time,int current_time);

  virtual double get_decay_rate(ImmuneSystemInformation h_immune_system_information,const int &age = 0) const = 0;

  virtual double get_acquire_rate(ImmuneSystemInformation h_immune_system_information,int latest_update_time,const int &age = 0) const = 0;

  virtual int type() const {
      return 0;
  }

};

#endif    /* IMMUNECOMPONENT_CUH */
