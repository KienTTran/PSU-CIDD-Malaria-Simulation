/* 
 * File:   ParasiteUpdateFunction.h
 * Author: Merlin
 *
 * Created on July 11, 2013, 2:51 PM
 */

#ifndef PARASITEDENSITYUPDATEFUNCTION_CUH
#define PARASITEDENSITYUPDATEFUNCTION_CUH

#include "Core/TypeDef.h"

namespace GPU{
    class ParasiteDensityUpdateFunction;
    class ClonalParasitePopulation;
    class Genotype;
}

class GPU::ParasiteDensityUpdateFunction {
 public:
    ParasiteDensityUpdateFunction();

  //    ParasiteUpdateFunction(const ParasiteUpdateFunction& orig);
  virtual ~ParasiteDensityUpdateFunction();

  virtual double get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) = 0;
  virtual int type() const {
      return 0;
    }

 private:

};

#endif    /* PARASITEDENSITYUPDATEFUNCTION_CUH */

