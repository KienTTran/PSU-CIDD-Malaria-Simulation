/* 
 * File:   PersonIndexByLocationMovingLevel.h
 * Author: Merlin
 *
 * Created on August 1, 2013, 9:26 PM
 */

#ifndef PERSONINDEXBYLOCATIONMOVINGLEVEL_CUH
#define    PERSONINDEXBYLOCATIONMOVINGLEVEL_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "PersonIndex.cuh"
#include "Gpu/Population/Person.cuh"

namespace GPU {
    class PersonIndexByLocationMovingLevel;
}

class GPU::PersonIndexByLocationMovingLevel : public GPU::PersonIndex {
 PROPERTY_REF(GPUPersonPtrVector3, vPerson);
 public:
  PersonIndexByLocationMovingLevel(const int &no_location = 1, const int &no_level = 1);

  //    PersonIndexByLocationMovingLevel(const PersonIndexByLocationMovingLevel& orig);
  virtual ~PersonIndexByLocationMovingLevel();

  void Initialize(const int &no_location = 1, const int &no_level = 1);

  virtual void add(GPU::Person *p);

  virtual void remove(GPU::Person *p);

  virtual std::size_t size() const;

  virtual void update();

  virtual void notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue);

 private:
  void remove_without_set_index(GPU::Person *p);

  void add(GPU::Person *p, const int &location, const int &moving_level);

  void change_property(GPU::Person *p, const int &location, const int &bitting_level);

};

#endif    /* PERSONINDEXBYLOCATIONMOVINGLEVEL_H */

