/* 
 * File:   GPU::PersonIndexByLocationStateAgeClass.h
 * Author: nguyentran
 *
 * Created on May 2, 2013, 10:41 AM
 */

#ifndef PERSONINDEXBYLOCATIONSTATEAGECLASS_CUH
#define    PERSONINDEXBYLOCATIONSTATEAGECLASS_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "PersonIndex.cuh"
#include "Gpu/Population/Person.cuh"

namespace GPU{
  class Person;
  class PersonIndexByLocationStateAgeClass;
};

class GPU::PersonIndexByLocationStateAgeClass : public GPU::PersonIndex {
 PROPERTY_REF(GPUPersonPtrVector4, vPerson);

 public:
  //    GPU::PersonIndexByLocationStateAgeClass();
  PersonIndexByLocationStateAgeClass(const int &no_location = 1, const int &no_host_state = 1,
                                     const int &no_age_class = 1);

  //    GPU::PersonIndexByLocationStateAgeClass(const GPU::PersonIndexByLocationStateAgeClass& orig);
  virtual ~PersonIndexByLocationStateAgeClass();

  void Initialize(const int &no_location = 1, const int &no_host_state = 1, const int &no_age_class = 1);

  virtual void add(GPU::Person *p);

  virtual void remove(GPU::Person *p);

  virtual std::size_t size() const;

  virtual void update();

  virtual void notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue);

 private:
  void remove_without_set_index(GPU::Person *p);

  void add(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);

  void change_property(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);
};

#endif    /* PERSONINDEXBYLOCATIONSTATEAGECLASS_H */

