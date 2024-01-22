/* 
 * File:   BloodParasite.h
 * Author: Merlin
 *
 * Created on July 11, 2013, 2:21 PM
 */

#ifndef CLONALPARASITEPOPULATION_CUH
#define CLONALPARASITEPOPULATION_CUH

#include "Core/PropertyMacro.h"
#include "Therapies/DrugType.h"
#include "Core/ObjectPool.h"
#include "Gpu/Population/Properties/IndexHandler.cuh"
#include "SingleHostClonalParasitePopulations.cuh"
#include "ClinicalUpdateFunction.cuh"
#include "ImmunityClearanceUpdateFunction.cuh"
#include "ParasiteDensityUpdateFunction.cuh"

class Therapy;

namespace GPU{
    class ClonalParasitePopulation;
    class SingleHostClonalParasitePopulations;
    class ClinicalUpdateFunction;
    class Genotype;
    class ImmuneSystem;
}

class GPU::ClonalParasitePopulation : public IndexHandler {
public:
  PROPERTY_HEADER(double, last_update_log10_parasite_density)

 PROPERTY_HEADER(double, gametocyte_level)
//    PROPERTY_REF(double, clearance_rate)
 PROPERTY_REF(int, first_date_in_blood)

 POINTER_PROPERTY(GPU::SingleHostClonalParasitePopulations, parasite_population)

 POINTER_PROPERTY_HEADER(GPU::Genotype, genotype)

 POINTER_PROPERTY(GPU::ParasiteDensityUpdateFunction, update_function)

 public:
    static const double LOG_ZERO_PARASITE_DENSITY;


 public:
    ClonalParasitePopulation(GPU::Genotype *genotype = nullptr);

  //    BloodParasite(const BloodParasite& orig);
  virtual ~ClonalParasitePopulation();

  __host__ double get_current_parasite_density(const int &current_time);

    __device__ __host__ double get_log10_infectious_density() const;

  bool resist_to(const int &drug_id) const;

  void update();

  void perform_drug_action(const double &percent_parasite_remove);

public:/* FOr GPU */
    __host__ void set_gpu_update_function(GPU::ParasiteDensityUpdateFunction *h_function);

    __device__ double get_current_parasite_density_gpu(GPU::ParasiteDensityUpdateFunction* d_update_function,
                                                       GPU::Genotype* d_genotype,
                                                       GPU::ImmuneSystem* d_immune_system,
                                                       ParasiteDensityLevel h_parasite_density_level,
                                                       ImmuneSystemInformation* d_immune_system_information,
                                                       const int current_time,
                                                       const int latest_update_time);
    double test_ = 2.0;
    __device__ __host__ double test(){
      test_ = 9.0;
      return test_;
    }

    __device__ __host__ double test2();

};

#endif    /* CLONALPARASITEPOPULATION_CUH */

