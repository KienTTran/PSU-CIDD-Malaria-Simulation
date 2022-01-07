/*
 * File:   IntParasiteDatabase.h
 * Author: Merlin
 *
 * Created on March 18, 2014, 3:06 PM
 */

#ifndef INTPARASITEDATABASE_H
#define INTPARASITEDATABASE_H

#include <map>

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "Genotype.h"

class Genotype;
class Config;

typedef std::map<ul, Genotype*> GenotypePtrMap;
typedef std::vector<std::vector<std::vector<double>>> MatingMatrix;

class GenotypeDatabase : public GenotypePtrMap {
  DISALLOW_COPY_AND_ASSIGN(GenotypeDatabase)

  DISALLOW_MOVE(GenotypeDatabase)

  VIRTUAL_PROPERTY_REF(MatingMatrix, mating_matrix)

  VIRTUAL_PROPERTY_REF(IntVector, weight)
public:
  std::map<std::string, unsigned int> aa_sequence_id_map;

public:
  GenotypeDatabase();

  virtual ~GenotypeDatabase();

  void add(Genotype* genotype);

  unsigned int get_id(const std::string& aa_sequence, Config* config);

  Genotype* get_genotype_from_alleles_structure(const IntVector& alleles);

private:
  unsigned int auto_id;
};

#endif /* INTPARASITEDATABASE_H */
