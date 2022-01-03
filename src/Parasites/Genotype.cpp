/*
 * File:   Genotype.cpp
 * Author: Merlin
 *
 * Created on March 17, 2014, 2:33 PM
 */
#include "Genotype.h"

#include <algorithm>

#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Model.h"
#include "Therapies/DrugDatabase.h"
#include "Therapies/SCTherapy.h"

Genotype::Genotype(const std::string &in_aa_sequence) : aa_sequence(in_aa_sequence) {
  // create aa structure
  std::string chromosome_str;
  std::istringstream tokenStream(in_aa_sequence);
  auto ii = 0;
  while (std::getline(tokenStream, chromosome_str, '|')) {
    std::string gene_str;
    std::istringstream chromosomeTokenStream(chromosome_str);
    auto jj = 0;
    while (std::getline(chromosomeTokenStream, gene_str, ',')) {
      aa_structure[ii].push_back(gene_str);
      jj++;
    }
    ii++;
  }

  // check if aa_sequence is valid

  // calculate cost of resistance

  // calculate ec50
}

Genotype::~Genotype() = default;

bool Genotype::resist_to(DrugType *dt) {
  // TODO: rework on this
  return false;
}

bool Genotype::resist_to(Therapy *therapy) {
  // TODO: rework on this
  return false;
}

Genotype *Genotype::combine_mutation_to(const int &locus, const int &value) {
  // TODO: rework on this
  return this;
}

double Genotype::get_EC50_power_n(DrugType *dt) const {
  return get_EC50(dt->id());
}

double Genotype::get_EC50(const int &drug_id) const {
  return Model::CONFIG->EC50_power_n_table()[genotype_id_][drug_id];
}

int Genotype::select_mutation_allele(const int &mutation_locus) {
  // TODO: rework on this
  return 0;
}

std::ostream &operator<<(std::ostream &os, const Genotype &e) {
  os << e.genotype_id_ << "\t";
  os << e.get_aa_sequence();
  return os;
}

std::string Genotype::get_aa_sequence() const {
  return aa_sequence;
  //  std::stringstream ss;
  //
  //  for (auto &chromosome : aa_structure) {
  //    ss << chromosome;
  //    if (&chromosome != &aa_structure.back()){
  //      ss << "|";
  //    }
  //  }
  //
  //  return ss.str();
}
bool Genotype::is_valid(const PfGeneInfo &gene_info) {
  for (int chromosome_i = 0; chromosome_i < 14; ++chromosome_i) {
    auto chromosome_info = gene_info.chromosome_infos[chromosome_i];
    if (chromosome_info.gene_infos.size() != aa_structure[chromosome_i].size()) return false;

    for (int gene_i = 0; gene_i < aa_structure[chromosome_i].size(); ++gene_i) {
      auto gene_info = chromosome_info.gene_infos[gene_i];
      auto max_aa_pos = gene_info.max_copy > 1 ? aa_structure[chromosome_i][gene_i].size() - 1
                                               : aa_structure[chromosome_i][gene_i].size();

      // check same size with aa postions info
      if (gene_info.aa_position_infos.size() != max_aa_pos) {
        std::cout << aa_structure[chromosome_i][gene_i] << std::endl;
        std::cout << gene_info.aa_position_infos.size() << std::endl;
        return false;
      }

      for (int aa_i = 0; aa_i < max_aa_pos; ++aa_i) {
        auto aa_pos_info = gene_info.aa_position_infos[aa_i];
        auto element = aa_structure[chromosome_i][gene_i][aa_i];

        if (std::find(aa_pos_info.amino_acids.begin(), aa_pos_info.amino_acids.end(), element)
            == aa_pos_info.amino_acids.end())
          return false;
      }

      // check number copy valid or not
      if (gene_info.max_copy > 1) {
        auto copy_number = (int)aa_structure[chromosome_i][gene_i].back() - 48;
        if (copy_number > gene_info.max_copy) {
          return false;
        }
      }
    }
  }

  return true;
}
