/* 
 * File:   ImmunityClearanceUpdateFunction.h
 * Author: Merlin
 *
 * Created on July 29, 2013, 5:49 PM
 */

#ifndef IMMUNITYCLEARANCEUPDATEFUNCTION_CUH
#define    IMMUNITYCLEARANCEUPDATEFUNCTION_CUH

#include "ParasiteDensityUpdateFunction.cuh"
#include "Core/PropertyMacro.h"

class Model;

namespace GPU{
    class ImmunityClearanceUpdateFunction;
    class ClonalParasitePopulation;
    class Genotype;
}

class GPU::ImmunityClearanceUpdateFunction : public GPU::ParasiteDensityUpdateFunction {
 POINTER_PROPERTY(Model, model)

 public:
  explicit __device__ __host__ ImmunityClearanceUpdateFunction(Model *model = nullptr);

  //    ImmunityClearanceUpdateFunction(const ImmunityClearanceUpdateFunction& orig);
  virtual ~ImmunityClearanceUpdateFunction();

  __host__ double get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) override;
  __device__ double get_current_parasite_density_gpu(ParasiteDensityLevel h_parasite_density_level,
                                                     double person_immune_parasite_size_after_t_days,
                                                     int duration) override;
    __device__ __host__  int type() const override {
    return 2;
  }

};

#endif    /* IMMUNITYCLEARANCEUPDATEFUNCTION_CUH */
