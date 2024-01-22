/*
 * File:   ParasitePopulation.h
 * Author: Merlin
 *
 * Created on July 11, 2013, 1:53 PM
 */

#ifndef SINGLEHOSTCLONALPARASITEPOPULATIONS_CUH
#define SINGLEHOSTCLONALPARASITEPOPULATIONS_CUH

#include <vector>

#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "ClonalParasitePopulation.cuh"

namespace GPU{
    class ClonalParasitePopulation;
    class SingleHostClonalParasitePopulations;
    class ParasiteDensityUpdateFunction;
    class Person;
    class DrugsInBlood;
}

class DrugType;

class GPU::SingleHostClonalParasitePopulations {

public:
  POINTER_PROPERTY(GPU::Person, person)

  POINTER_PROPERTY(std::vector<GPU::ClonalParasitePopulation *>, parasites)

public:
  // this value will be automatically updated daily in the function clear_cured_parasites
  // in order to have accurate sum of all density
  double log10_total_infectious_denstiy{-1000.0};

public:
    __device__ __host__ SingleHostClonalParasitePopulations(GPU::Person *person = nullptr);

  //    ParasitePopulation(const ParasitePopulation& orig);
  virtual ~SingleHostClonalParasitePopulations();

  void init();

  virtual int size();

  virtual void add(GPU::ClonalParasitePopulation *blood_parasite);

  virtual void remove(GPU::ClonalParasitePopulation *blood_parasite);

  virtual void remove(const int &index);

  virtual int latest_update_time() const;

  virtual bool contain(GPU::ClonalParasitePopulation *blood_parasite);

  void change_all_parasite_update_function(GPU::ParasiteDensityUpdateFunction *from,
                                           GPU::ParasiteDensityUpdateFunction *to) const;

  void update();

  void clear_cured_parasites();

  void clear();

  void update_by_drugs(GPU::DrugsInBlood *drugs_in_blood) const;

  bool has_detectable_parasite() const;

  bool is_gametocytaemic() const;
};

#endif /* SINGLEHOSTCLONALPARASITEPOPULATIONS_H */
