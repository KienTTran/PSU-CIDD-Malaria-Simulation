/* 
 * File:   PersonIndex.h
 * Author: nguyentran
 *
 * Created on April 17, 2013, 10:01 AM
 */

#ifndef PERSONINDEX_CUH
#define PERSONINDEX_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Population/Person.cuh"

namespace GPU{
    class PersonIndex;
    class Person;
}

class GPU::PersonIndex {
 public:
  PersonIndex();

  virtual ~PersonIndex();

  virtual void add(GPU::Person *p) = 0;

  virtual void remove(GPU::Person *p) = 0;

  virtual std::size_t size() const = 0;

  virtual void update() = 0;

  virtual void
  notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue) = 0;

 private:

};

#endif    /* PERSONINDEX_H */

