/*
 * File:   DrugType.h
 * Author: nguyentran
 *
 * Created on June 3, 2013, 10:59 AM
 */

#ifndef DRUGTYPE_CUH
#define DRUGTYPE_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

namespace GPU{
    class Genotype;
    class DrugType;
}

class Config;

class GPU::DrugType {
public:
    struct ResistantAALocation {
        int chromosome_id { -1 };
        int gene_id { -1 };
        int aa_id { -1 };
        int aa_index_in_aa_string { -1 };
        bool is_copy_number { false };
    };

public:
  DISALLOW_COPY_AND_ASSIGN(DrugType)

  VIRTUAL_PROPERTY_REF(int, id)

  VIRTUAL_PROPERTY_REF(std::string, name)

  VIRTUAL_PROPERTY_REF(double, drug_half_life)

  VIRTUAL_PROPERTY_REF(double, maximum_parasite_killing_rate)

  VIRTUAL_PROPERTY_REF(DoubleVector, age_group_specific_drug_concentration_sd);

  VIRTUAL_PROPERTY_REF(DoubleVector, age_specific_drug_absorption);

  VIRTUAL_PROPERTY_REF(double, k)

  VIRTUAL_PROPERTY_REF(double, cut_off_percent)

public:
  double base_EC50 { 0 };
  std::vector<ResistantAALocation> resistant_aa_locations {};

public:
  DrugType();

  virtual ~DrugType();

  virtual double get_parasite_killing_rate_by_concentration(const double &concentration, const double &EC50_power_n);

  virtual double n();

  virtual void set_n(const double &n);

  int get_total_duration_of_drug_activity(const int &dosing_days) const;

  void populate_resistant_aa_locations(Config *config_);

private:
  double n_;
  //    double EC50_;
};

#endif /* DRUGTYPE_CUH */
