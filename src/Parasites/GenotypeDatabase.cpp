/*
 * File:   IntParasiteDatabase.cpp
 * Author: Merlin
 *
 * Created on March 18, 2014, 3:06 PM
 */

#include "GenotypeDatabase.h"

#include "Core/Config/Config.h"
#include "Genotype.h"
#include "Helpers/NumberHelpers.h"
#include "Model.h"
#include "Core/TypeDef.h"

GenotypeDatabase::GenotypeDatabase() : auto_id(0) {}

GenotypeDatabase::~GenotypeDatabase() {
  for (auto &i : *this) {
    delete i.second;
  }
  clear();
}

void GenotypeDatabase::add(Genotype *genotype) {
  if (this->find(genotype->genotype_id) != this->end()) {
    delete (*this)[genotype->genotype_id];
  }
  (*this)[genotype->genotype_id] = genotype;
}

Genotype *GenotypeDatabase::get_genotype_from_alleles_structure(const IntVector &alleles) {
  auto id = 0;
  for (auto i = 0; i < alleles.size(); i++) {
    id += weight_[i] * alleles[i];
  }
  return this->at(id);
}

unsigned int GenotypeDatabase::get_id(const std::string &aa_sequence, Config *config) {
  return get_genotype(aa_sequence, config)->genotype_id;
}
Genotype *GenotypeDatabase::get_genotype(const std::string &aa_sequence, Config *config) {
  if (aa_sequence_id_map.find(aa_sequence) == aa_sequence_id_map.end()) {
    // not yet exist then initialize new genotype
    auto new_id = auto_id;
    auto new_genotype = new Genotype(aa_sequence);
    new_genotype->genotype_id = new_id;
    new_genotype->resistant_recombinations_in_mosquito = std::vector<MosquitoRecombinedGenotypeInfo>();

    // check if aa_sequence is valid
    if (!new_genotype->is_valid(config->pf_genotype_info())) {
      LOG(FATAL) << "Invalid genotype: " << aa_sequence;
    }

    // calculate cost of resistance
    new_genotype->calculate_daily_fitness(config->pf_genotype_info());

    // calculate ec50
    new_genotype->calculate_EC50_power_n(config->pf_genotype_info(), config->drug_db());

    new_genotype->override_EC50_power_n(config->override_ec50_patterns(), config->drug_db());

    // add min ec50 of each drug to db
    for(int drug_id = 0; drug_id <  config->drug_db()->size(); drug_id++){
      if (drug_id_ec50.find(drug_id) == drug_id_ec50.end()){
        if (drug_id_ec50[drug_id].find(new_genotype->get_aa_sequence()) == drug_id_ec50[drug_id].end()){
          drug_id_ec50[drug_id][new_genotype->get_aa_sequence()] = new_genotype->get_EC50_power_n(config->drug_db()->at(drug_id));
        }
        else{
          if (drug_id_ec50[drug_id][new_genotype->get_aa_sequence()] > new_genotype->get_EC50_power_n(config->drug_db()->at(drug_id))){
            drug_id_ec50[drug_id][new_genotype->get_aa_sequence()] = new_genotype->get_EC50_power_n(config->drug_db()->at(drug_id));
          }
        }
      }
    }

    aa_sequence_id_map[aa_sequence] = new_genotype;
    add(new_genotype);
    new_genotype->resistant_recombinations_in_mosquito = std::vector<MosquitoRecombinedGenotypeInfo>();

    auto_id++;
    return new_genotype;
  } else {
    return aa_sequence_id_map[aa_sequence];
  }
}

double GenotypeDatabase::get_min_ec50(int drug_id){
  auto it = min_element(drug_id_ec50[drug_id].begin(), drug_id_ec50[drug_id].end(),
                                                 [](const auto& l, const auto& r) { return l.second < r.second; });

  return it->second;
}

