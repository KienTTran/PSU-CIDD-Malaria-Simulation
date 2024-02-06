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
  ImmunityClearanceUpdateFunction(Model *model = nullptr);

  //    ImmunityClearanceUpdateFunction(const ImmunityClearanceUpdateFunction& orig);
  virtual ~ImmunityClearanceUpdateFunction();

    double get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) override;
     int type() const override {
    return 2;
  }

};

#endif    /* IMMUNITYCLEARANCEUPDATEFUNCTION_CUH */
