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
    /*
     * https://stackoverflow.com/questions/307765/how-do-i-check-if-an-objects-type-is-a-particular-subclass-in-c
     * For checking type of class to implement polymorphism in CUDA
    */
    __device__ __host__ ParasiteDensityUpdateFunction();

  //    ParasiteUpdateFunction(const ParasiteUpdateFunction& orig);
  virtual ~ParasiteDensityUpdateFunction();

    __host__ virtual double get_current_parasite_density(GPU::ClonalParasitePopulation *parasite, int duration) = 0;
    __device__ virtual double get_current_parasite_density_gpu(ParasiteDensityLevel h_parasite_density_level,
                                                               double person_immune_parasite_size_after_t_days,
                                                               int duration) = 0;
    __device__ __host__ virtual int type() const {
      return 0;
    }

 private:

};

#endif    /* PARASITEDENSITYUPDATEFUNCTION_CUH */

