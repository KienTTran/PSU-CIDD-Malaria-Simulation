/*
 * File:   IntParasiteDatabase.h
 * Author: Merlin
 *
 * Created on March 18, 2014, 3:06 PM
 */

#ifndef GENOTYPEDATABASE_CUH
#define GENOTYPEDATABASE_CUH

#include <map>

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "Genotype.cuh"

namespace GPU{
    class GenotypeDatabase;
    class Genotype;
    typedef std::map<ul, GPU::Genotype*> GenotypePtrMap;
}
class Config;


class GPU::GenotypeDatabase : public GPU::GenotypePtrMap {
//  VIRTUAL_PROPERTY_REF(MatingMatrix, mating_matrix)

  VIRTUAL_PROPERTY_REF(IntVector, weight)

public:
  std::map<std::string, GPU::Genotype*> aa_sequence_id_map;

public:
  GenotypeDatabase();

  virtual ~GenotypeDatabase();

  void add(GPU::Genotype* genotype);

  GPU::Genotype* get_genotype(const std::string& aa_sequence, Config* pConfig);

  unsigned int get_id(const std::string& aa_sequence, Config* config);

  GPU::Genotype* get_genotype_from_alleles_structure(const IntVector& alleles);

private:
  unsigned int auto_id {0};
};

#endif /* GENOTYPEDATABASE_CUH */
