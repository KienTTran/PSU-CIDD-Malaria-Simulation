//
// Created by nguyentd on 3/11/2022.
//

#ifndef POMS_SRC_MOSQUITO_MOSQUITO_H
#define POMS_SRC_MOSQUITO_MOSQUITO_H
#include <vector>

#include "Core/Config/Config.h"
#include "Core/PropertyMacro.h"
#include "Therapies/SCTherapy.h"

class Genotype;
class Model;
class Config;
class Population;

class Mosquito {
public:
  Mosquito(const Mosquito &) = delete;
  const Mosquito &operator=(const Mosquito &) = delete;

public:
  Model *model { nullptr };

public:
    std::map<int,double> drug_id_min_ec50;
    std::vector<std::vector<std::string>> double_resistant_list = {{"AL:2-2", "AL:2-3", "AL:2-4"},
                                                                   {"ASAQ:2-2", "ASAQ:2-3", "ASAQ:2-4"},
                                                                   {"DHA-PPQ:2-2"}};

public:
  explicit Mosquito(Model *model = nullptr);
  virtual ~Mosquito() = default;

  void initialize(Config *config);

  void infect_new_cohort_in_PRMC(Config *config, Random *random, Population *population, const int &tracking_index);

public:
  std::vector<std::vector<std::vector<Genotype *>>> genotypes_table; /* Mosquito table */

  [[nodiscard]] static std::vector<unsigned int> build_interrupted_feeding_indices(
      Random *random, const double &interrupted_feeding_rate, const int &prmc_size);

  int random_genotype(int location, int tracking_index);

  // this function will populate values for both parasite densities and genotypes that carried by a person
  void get_genotypes_profile_from_person(Person *person, std::vector<Genotype *> &sampling_genotypes,
                                         std::vector<double> &relative_infectivity_each_pp);

  bool genotype_resistant_to(Genotype *genotype, std::string resistance, int therapy_id);

  bool genotype_resistant_to(Config* config, std::vector<Genotype*> parent_genotypes, SCTherapy * sc_therapy, Genotype *genotype, std::string resistance, int therapy_id, bool verbose = false);

  std::vector<std::string> split_string(std::string str, char delimiter);
};

#endif  // POMS_SRC_MOSQUITO_MOSQUITO_H
