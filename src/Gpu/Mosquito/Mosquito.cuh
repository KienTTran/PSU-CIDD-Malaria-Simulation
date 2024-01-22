//
// Created by nguyentd on 3/11/2022.
//

#ifndef POMS_SRC_MOSQUITO_MOSQUITO_CUH
#define POMS_SRC_MOSQUITO_MOSQUITO_CUH
#include <vector>

#include "Core/Config/Config.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class Genotype;
    class Person;
    class Mosquito;
    class Population;
}
class Config;
class Random;

class GPU::Mosquito {
public:
  Mosquito(const Mosquito &) = delete;
  const Mosquito &operator=(const Mosquito &) = delete;

public:
  Model *model { nullptr };

public:
  explicit Mosquito(Model *model = nullptr);
  virtual ~Mosquito() = default;

  void initialize(Config *config);

  void infect_new_cohort_in_PRMC(Config *config, ::Random *random, GPU::Population *population, const int &tracking_index);

public:
  std::vector<std::vector<std::vector<GPU::Genotype *>>> genotypes_table; /* Mosquito table */

  [[nodiscard]] static std::vector<unsigned int> build_interrupted_feeding_indices(
      ::Random *random, const double &interrupted_feeding_rate, const int &prmc_size);

  int random_genotype(int location, int tracking_index);

  // this function will populate values for both parasite densities and genotypes that carried by a person
  void get_genotypes_profile_from_person(GPU::Person *person, std::vector<GPU::Genotype *> &sampling_genotypes,
                                         std::vector<double> &relative_infectivity_each_pp);
};

#endif  // POMS_SRC_MOSQUITO_MOSQUITO_H
