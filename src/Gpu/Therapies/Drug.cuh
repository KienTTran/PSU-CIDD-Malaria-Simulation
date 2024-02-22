/*
 * File:   Drug.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 1:47 PM
 */

#ifndef DRUG_CUH
#define DRUG_CUH

#include "Core/PropertyMacro.h"
#include "Core/ObjectPool.h"
#include "DrugType.cuh"

namespace GPU {
  class DrugsInBlood;
  class Drug;
}

class GPU::Drug {

 DISALLOW_COPY_AND_ASSIGN(Drug)

 PROPERTY_REF(int, dosing_days)

 PROPERTY_REF(int, start_time)

 PROPERTY_REF(int, end_time)

 PROPERTY_REF(double, last_update_value)

 PROPERTY_REF(int, last_update_time)

 PROPERTY_REF(double, starting_value)

 POINTER_PROPERTY(GPU::DrugType, drug_type)

 POINTER_PROPERTY(GPU::DrugsInBlood, person_drugs)

 public:
  explicit Drug(GPU::DrugType *drug_type = nullptr);

  //    Drug(const Drug& orig);
  virtual ~Drug();

  void update();

  double get_current_drug_concentration(int currentTime);

  double get_mutation_probability() const;

  double get_mutation_probability(double currentDrugConcentration) const;

  void set_number_of_dosing_days(int dosingDays);

  double get_parasite_killing_rate(const int &genotype_id) const;
};

#endif    /* DRUG_CUH */