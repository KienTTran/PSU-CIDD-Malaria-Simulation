//
// Created by nguyentd on 3/11/2022.
//

#include "Mosquito.h"

#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Model.h"

Mosquito::Mosquito(Model *model) : model { model } {}

void Mosquito::initialize(Config *config) {
  genotypes_table.clear();

  genotypes_table = std::vector<std::vector<std::vector<Genotype *>>>(
      config->number_of_locations(),
      std::vector<std::vector<Genotype *>>(config->number_of_tracking_days(),
                                           std::vector<Genotype *>(config->mosquito_config().prmc_size, nullptr))

  );
  VLOG(0) << "[Mosquito] Size: " << config->mosquito_config().prmc_size;
  if (!config->mosquito_config().interrupted_feeding_rate_raster.empty()) {
    // read from raster
    //    MosquitoData::get_instance().load_raster_from_path(Model::CONFIG->mosquito_config().interrupted_feeding_rate_raster,
    //    MosquitoData::InteruptedFeedingRate);
    LOG(FATAL) << "Raster is not supported in version 3.x!!!";
  } else {
    if (config->mosquito_config().interrupted_feeding_rate.size() == 1) {
      double if_rate = config->mosquito_config().interrupted_feeding_rate[0];
      config->mosquito_config().interrupted_feeding_rate = std::vector<double>(config->number_of_locations(), if_rate);
    } else if (config->mosquito_config().interrupted_feeding_rate.size() != config->number_of_locations()) {
      LOG(FATAL) << "Number of element of interrupted feeding rate should be 1 or equal to number of locations!!!!";
    }
  }
}

void Mosquito::infect_new_cohort_in_PRMC(Config *config, Random *random, const int &tracking_index) {
  // for each location fill prmc at tracking_index row with sampling genotypes
  for (int location = 0; location < config->number_of_locations(); location++) {

    // multinomial sampling based on relative infectivity

    std::vector<unsigned int> all_interrupted_feeding = build_interrupted_feeding_indices(
        random, config->mosquito_config().interrupted_feeding_rate[location], config->mosquito_config().prmc_size);

    // uniform sampling in all person


    // recombination
    // *p1 , *p2, bool is_interrupted  ===> *genotype

  }
}

std::vector<unsigned int> Mosquito::build_interrupted_feeding_indices(Random *random,
                                                                      const double &interrupted_feeding_rate,
                                                                      const int &prmc_size) {
  int number_of_interrupted_feeding = random->random_poisson(interrupted_feeding_rate * prmc_size);

  std::vector<unsigned int> all_interrupted_feeding(number_of_interrupted_feeding, 1);
  all_interrupted_feeding.resize(prmc_size, 0);

  random->random_shuffle(&all_interrupted_feeding[0], all_interrupted_feeding.size(), sizeof(unsigned int));
  return all_interrupted_feeding;
}
