/*
 * File:   IntParasiteDatabase.cpp
 * Author: Merlin
 *
 * Created on March 18, 2014, 3:06 PM
 */

#include "GenotypeDatabase.cuh"

#include "Core/Config/Config.h"
#include "Genotype.cuh"
#include "Helpers/NumberHelpers.h"
#include "Model.h"

GPU::GenotypeDatabase::GenotypeDatabase() : auto_id(0) {}

GPU::GenotypeDatabase::~GenotypeDatabase() {
  for (auto &i : *this) {
    delete i.second;
  }
  clear();
}

void GPU::GenotypeDatabase::add(GPU::Genotype *genotype) {
  if (this->find(genotype->genotype_id) != this->end()) {
    delete (*this)[genotype->genotype_id];
  }
  (*this)[genotype->genotype_id] = genotype;
}

GPU::Genotype *GPU::GenotypeDatabase::get_genotype_from_alleles_structure(const IntVector &alleles) {
  auto id = 0;
  for (auto i = 0; i < alleles.size(); i++) {
    id += weight_[i] * alleles[i];
  }
  return this->at(id);
}

unsigned int GPU::GenotypeDatabase::get_id(const std::string &aa_sequence, Config *config) {
  return get_genotype(aa_sequence, config)->genotype_id;
}
GPU::Genotype *GPU::GenotypeDatabase::get_genotype(const std::string &aa_sequence, Config *config) {
  if (aa_sequence_id_map.find(aa_sequence) == aa_sequence_id_map.end()) {
    // not yet exist then initialize new genotype
    auto new_id = auto_id;
    auto new_genotype = new GPU::Genotype(aa_sequence);
    new_genotype->genotype_id = new_id;

    // check if aa_sequence is valid
    if (!new_genotype->is_valid(config->pf_genotype_info())) {
      LOG(FATAL) << "Invalid genotype: " << aa_sequence;
    }

    // calculate cost of resistance
    new_genotype->calculate_daily_fitness(config->pf_genotype_info());

    // calculate ec50
    new_genotype->calculate_EC50_power_n(config->pf_genotype_info(), config->drug_db());

    new_genotype->override_EC50_power_n(config->override_ec50_patterns(), config->drug_db());

    aa_sequence_id_map[aa_sequence] = new_genotype;
    add(new_genotype);

    auto_id++;
    return new_genotype;
  } else {
    return aa_sequence_id_map[aa_sequence];
  }
}
