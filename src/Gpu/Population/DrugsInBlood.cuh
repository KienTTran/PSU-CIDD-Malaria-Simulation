/* 
 * File:   DrugsInBlood.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 1:47 PM
 */

#ifndef DRUGSINBLOOD_CUH
#define    DRUGSINBLOOD_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "Therapies/DrugType.h"

namespace GPU{
    class Person;
    class DrugsInBlood;
}

class GPU::DrugsInBlood {
 POINTER_PROPERTY(GPU::Person, person)

 POINTER_PROPERTY(DrugPtrMap, drugs)

 public:
  explicit DrugsInBlood(GPU::Person *person = nullptr);

  //    DrugsInBlood(const DrugsInBlood& orig);
  virtual ~DrugsInBlood();

  void init();

  Drug *add_drug(Drug *drug);

  bool is_drug_in_blood(DrugType *drug_type) const;

  bool is_drug_in_blood(int drug_type_id) const;

  void remove_drug(Drug *drug) const;

  void remove_drug(const int &drug_type_id) const;

  Drug *get_drug(const int &type_id) const;

  std::size_t size() const;

  void clear() const;

  void update() const;

  void clear_cut_off_drugs_by_event(GPU::Event *event) const;

};

#endif    /* DRUGSINBLOOD_H */
