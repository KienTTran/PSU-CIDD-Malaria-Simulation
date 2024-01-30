/* 
 * File:   BloodParasite.h
 * Author: Merlin
 *
 * Created on July 11, 2013, 2:21 PM
 */

#ifndef CLONALPARASITEPOPULATION_CUH
#define CLONALPARASITEPOPULATION_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Therapies/DrugType.cuh"
#include "Core/ObjectPool.h"
#include "Gpu/Population/Properties/IndexHandler.cuh"
#include "SingleHostClonalParasitePopulations.cuh"
#include "ClinicalUpdateFunction.cuh"
#include "ImmunityClearanceUpdateFunction.cuh"
#include "ParasiteDensityUpdateFunction.cuh"


namespace GPU{
    class ClonalParasitePopulation;
    class SingleHostClonalParasitePopulations;
    class ClinicalUpdateFunction;
    class Genotype;
    class ImmuneSystem;
    class Therapy;
    class Person;
}

class GPU::ClonalParasitePopulation : public GPU::IndexHandler {
public:
  PROPERTY_HEADER(double, last_update_log10_parasite_density)

 PROPERTY_HEADER(double, gametocyte_level)
//    PROPERTY_REF(double, clearance_rate)
 PROPERTY_REF(int, first_date_in_blood)

 POINTER_PROPERTY(GPU::SingleHostClonalParasitePopulations, parasite_population)

 POINTER_PROPERTY(GPU::Person, person)

 POINTER_PROPERTY_HEADER(GPU::Genotype, genotype)

 POINTER_PROPERTY_HEADER(GPU::ParasiteDensityUpdateFunction, update_function)

 public:
    constexpr static const double LOG_ZERO_PARASITE_DENSITY = -1000.0;
    PROPERTY_REF(long, id);
    PROPERTY_REF(int, index);

 public:
    ClonalParasitePopulation(GPU::Genotype *genotype = nullptr);

  //    BloodParasite(const BloodParasite& orig);
  virtual ~ClonalParasitePopulation();

  __host__ double get_current_parasite_density(const int &current_time);

    __device__ __host__ double get_log10_infectious_density() const;

  bool resist_to(const int &drug_id) const;

  void update();

  void perform_drug_action(const double &percent_parasite_remove);


};

#endif    /* CLONALPARASITEPOPULATION_CUH */

