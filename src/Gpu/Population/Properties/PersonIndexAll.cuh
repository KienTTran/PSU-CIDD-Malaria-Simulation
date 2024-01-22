/* 
 * File:   PersonIndexAll.h
 * Author: nguyentran
 *
 * Created on April 17, 2013, 10:15 AM
 */

#ifndef PERSONINDEXALL_CUH
#define    PERSONINDEXALL_CUH

#include "PersonIndex.cuh"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

namespace GPU {
    class PersonIndexAll;
}

class GPU::PersonIndexAll : public GPU::PersonIndex {
 PROPERTY_REF(GPUPersonPtrVector, vPerson)

 public:
  PersonIndexAll();

  virtual ~PersonIndexAll();

  virtual void add(GPU::Person *p);

  virtual void remove(GPU::Person *p);

  virtual std::size_t size() const;

  virtual void update();

  virtual void notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue);

 private:

};

#endif    /* PERSONINDEXALL_H */

