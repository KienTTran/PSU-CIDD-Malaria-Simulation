//
// Created by nguyentd on 3/11/2022.
//

#include "Mosquito.h"

#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Core/TypeDef.h"
#include "Model.h"
#include "MDC/ModelDataCollector.h"
#include "Population/Population.h"
#include "Population/SingleHostClonalParasitePopulations.h"
#include "Therapies/SCTherapy.h"
#include <algorithm>
#include "easylogging++.h"

Mosquito::Mosquito(Model *model) : model { model } {}

void Mosquito::initialize(Config *config) {
  genotypes_table.clear();

  genotypes_table = std::vector<std::vector<std::vector<Genotype *>>>(
      config->number_of_tracking_days(),
      std::vector<std::vector<Genotype *>>(config->number_of_locations(),
                                           std::vector<Genotype *>(config->mosquito_config().prmc_size, nullptr)));
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
    // add min ec50 of each drug to db
    for (int drug_id = 0; drug_id < config->drug_db()->size(); drug_id++) {
        drug_id_min_ec50[drug_id] = pow(config->drug_db()->at(drug_id)->base_EC50, config->drug_db()->at(drug_id)->n());
        LOG(INFO) << fmt::format("Drug id: {} base_EC50: {} n: {} min_EC50: {}",drug_id,
               config->drug_db()->at(drug_id)->base_EC50, config->drug_db()->at(drug_id)->n(),
               drug_id_min_ec50[drug_id]);
    }
}

void Mosquito::infect_new_cohort_in_PRMC(Config *config, Random *random, Population *population,
                                         const int &tracking_index) {
  // for each location fill prmc at tracking_index row with sampling genotypes
  for (int loc = 0; loc < config->number_of_locations(); loc++) {
    LOG(TRACE) << "Day " << Model::SCHEDULER->current_time()
               << " ifr = " << Model::CONFIG->mosquito_config().interrupted_feeding_rate[loc];
    // if there is no parasites in location
    if (population->current_force_of_infection_by_location[loc] <= 0) {
      for (int i = 0; i < config->mosquito_config().prmc_size; ++i) {
        genotypes_table[tracking_index][loc][i] = nullptr;
      }
      return;
    }
    // multinomial sampling of people based on their relative infectivity (summing across all clones inside that person)
    auto first_sampling = random->roulette_sampling<Person>(config->mosquito_config().prmc_size,
                                                            population->individual_foi_by_location[loc],
                                                            population->all_alive_persons_by_location[loc], false,
                                                            population->current_force_of_infection_by_location[loc]);

    std::vector<unsigned int> interrupted_feeding_indices = build_interrupted_feeding_indices(
        random, config->mosquito_config().interrupted_feeding_rate[loc], config->mosquito_config().prmc_size);

    // uniform sampling in all person
    auto second_sampling = random->roulette_sampling<Person>(config->mosquito_config().prmc_size,
                                                             population->individual_relative_biting_by_location[loc],
                                                             population->all_alive_persons_by_location[loc], true);

    // recombination
    // *p1 , *p2, bool is_interrupted  ===> *genotype
    std::vector<Genotype *> sampled_genotypes;
    std::vector<double> relative_infectivity_each_pp;

    for (int if_index = 0; if_index < interrupted_feeding_indices.size(); ++if_index) {
      // clear() is used to avoid memory reallocation
      sampled_genotypes.clear();
      relative_infectivity_each_pp.clear();

      if (config->within_host_induced_free_recombination()) {
        // get all infectious parasites from first person
        get_genotypes_profile_from_person(first_sampling[if_index], sampled_genotypes, relative_infectivity_each_pp);

        if (sampled_genotypes.empty()) {
          LOG(FATAL) << "first person has no infectious parasites, log10_total_infectious_denstiy = "
                     << first_sampling[if_index]->all_clonal_parasite_populations()->log10_total_infectious_denstiy;
        }

        if (interrupted_feeding_indices[if_index]) {
          // if second person is the same as first person, re-select second person until it is different from first.
          // this is to avoid recombination between the same person because in this case the interrupted feeding is true,
          // this is worst case scenario
          auto temp_if = if_index;
          while (second_sampling[temp_if] == first_sampling[if_index]) {
            temp_if = random->random_uniform(second_sampling.size());
          }
          // interrupted feeding occurs
          get_genotypes_profile_from_person(second_sampling[temp_if], sampled_genotypes, relative_infectivity_each_pp);
          //Count interrupted feeding events with within host induced recombination on
          Model::DATA_COLLECTOR->mosquito_recombination_events_count()[loc][0]++;
        }

        if (sampled_genotypes.empty()) {
          LOG(FATAL) << "sampling_genotypes should not be empty";
        }
      } else {
        sampled_genotypes.clear();
        relative_infectivity_each_pp.clear();
        get_genotypes_profile_from_person(first_sampling[if_index], sampled_genotypes, relative_infectivity_each_pp);
        // get exactly 1 infectious parasite from first person
        auto first_genotype =
            random->roulette_sampling_tuple<Genotype>(1, relative_infectivity_each_pp, sampled_genotypes, false)[0];

        std::tuple<Genotype *, double> second_genotype = std::make_tuple(nullptr, 0.0);

        if (interrupted_feeding_indices[if_index]) {
          // if second person is the same as first person, re-select second person until it is different from first.
          // this is to avoid recombination between the same person because in this case the interrupted feeding is true,
          // this is worst case scenario
          auto temp_if = if_index;
          while (second_sampling[temp_if] == first_sampling[if_index]) {
            temp_if = random->random_uniform(second_sampling.size());
          }
          sampled_genotypes.clear();
          relative_infectivity_each_pp.clear();
          get_genotypes_profile_from_person(second_sampling[temp_if], sampled_genotypes, relative_infectivity_each_pp);

          if (sampled_genotypes.size() > 0) {
            second_genotype = random->roulette_sampling_tuple<Genotype>(1, relative_infectivity_each_pp,
                                                                        sampled_genotypes, false)[0];
          }
          //Count interrupted feeding events with within host induced recombination off
          Model::DATA_COLLECTOR->mosquito_recombination_events_count()[loc][0]++;
        }

        sampled_genotypes.clear();
        relative_infectivity_each_pp.clear();
        sampled_genotypes.push_back(std::get<0>(first_genotype));
        relative_infectivity_each_pp.push_back(std::get<1>(first_genotype));

        if (std::get<0>(second_genotype) != nullptr) {
          sampled_genotypes.push_back(std::get<0>(second_genotype));
          relative_infectivity_each_pp.push_back(std::get<1>(second_genotype));
        }
      }

      auto parent_genotypes = random->roulette_sampling<Genotype>(2, relative_infectivity_each_pp, sampled_genotypes, false);

      Genotype *sampled_genotype =
          (parent_genotypes[0]->aa_sequence == parent_genotypes[1]->aa_sequence)
              ? parent_genotypes[0]
              : Genotype::free_recombine(config, random, parent_genotypes[0], parent_genotypes[1]);

      genotypes_table[tracking_index][loc][if_index] = sampled_genotype;


      //Count DHA-PPQ(8) ASAQ(7) AL(6)
      //Count if male genotype resists to one drug and female genotype resists to another drug only, right now work on double and triple resistant only
      //when genotype ec50_power_n == min_ec50, it is sensitive to that drug
        if (Model::SCHEDULER->current_time() >= Model::CONFIG->start_of_comparison_period())
        {
            /*
             * Print our recombination for counting later
             * */
            auto resistant_tracker_info = std::make_tuple(Model::SCHEDULER->current_time(),parent_genotypes[0]->genotype_id, parent_genotypes[1]->genotype_id, sampled_genotype->genotype_id);
            if (std::find(Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_tracker[loc].begin(),
                          Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_tracker[loc].end(), resistant_tracker_info)
                == Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_tracker[loc].end()){
                Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_tracker[loc].push_back(resistant_tracker_info);
            }
        }
      //Count number of bites
      Model::DATA_COLLECTOR->mosquito_recombination_events_count()[loc][1]++;
    }
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

int Mosquito::random_genotype(int location, int tracking_index) {
    auto genotype_index = Model::RANDOM->random_uniform_int(0, Model::CONFIG->mosquito_config().prmc_size);
    return genotypes_table[tracking_index][location][genotype_index]->genotype_id;
}

void Mosquito::get_genotypes_profile_from_person(Person *person, std::vector<Genotype *> &sampling_genotypes,
                                                 std::vector<double> &relative_infectivity_each_pp) {
  for (auto *pp : *person->all_clonal_parasite_populations()->parasites()) {
    //Select parasites based on gametocyte density
    auto clonal_foi = pp->gametocyte_level() * Person::relative_infectivity(pp->last_update_log10_parasite_density());
    if (clonal_foi > 0) {
      relative_infectivity_each_pp.push_back(clonal_foi);
      sampling_genotypes.push_back(pp->genotype());
    }
  }
}

std::vector<std::string> Mosquito::split_string(std::string str, char delimiter) {
  std::vector<std::string> internal;
  std::stringstream ss(str); // Turn the string into a stream.
  std::string tok;
  while (getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
  return internal;
}

std::string Mosquito::get_old_genotype_string(std::string new_genotype){
    std::vector<std::string> pattern_chr = split_string(new_genotype,'|');
    std::string old_chr_5 = pattern_chr[6].substr(0, 1);
    std::string old_chr_7 = "";
    if(pattern_chr[4].substr(2, 1) == "2")
        old_chr_7 = pattern_chr[4].substr(0, 2)+pattern_chr[4].substr(0, 2);
    else
        old_chr_7 = pattern_chr[4].substr(0, 2)+"--";
    std::string old_chr_13 = pattern_chr[12].substr(10, 1);
    std::string old_chr_14 = pattern_chr[13].substr(0, 1);
    std::string old_chr_x = pattern_chr[6].substr(6, 1);
    return old_chr_5+old_chr_7+old_chr_13+old_chr_14;
}

