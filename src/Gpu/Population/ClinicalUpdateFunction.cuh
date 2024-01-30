/* 
 * File:   ClinicalUpdateFunction.h
 * Author: Merlin
 *
 * Created on July 29, 2013, 5:43 PM
 */

#ifndef CLINICALUPDATEFUNCTION_CUH
#define CLINICALUPDATEFUNCTION_CUH

#include "ParasiteDensityUpdateFunction.cuh"
#include "Core/PropertyMacro.h"

class Model;

namespace GPU{
    class ClinicalUpdateFunction;
    class Genotype;
}

class GPU::ClinicalUpdateFunction : public GPU::ParasiteDensityUpdateFunction {
 POINTER_PROPERTY(Model, model)

 public:
  ClinicalUpdateFunction(Model *model = nullptr);
  //    ClinicalUpdateFunction(const ClinicalUpdateFunction& orig);

  virtual ~ClinicalUpdateFunction();

    __device__ __host__ double get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) override;
    __device__ __host__ int type() const override {
        return 1;
  }


};

#endif    /* CLINICALUPDATEFUNCTION_CUH */
