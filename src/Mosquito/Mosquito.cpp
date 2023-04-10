//
// Created by nguyentd on 3/11/2022.
//

#include "Mosquito.h"

#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Model.h"
#include "MDC/ModelDataCollector.h"
#include "Population/Population.h"
#include "Population/SingleHostClonalParasitePopulations.h"
#include "Therapies/SCTherapy.h"
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
        printf("Drug id: %d base_EC50: %.3f n: %.1f min_EC50: %.10f\n", drug_id,
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
    // multinomial sampling based on relative infectivity
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
    std::vector<Genotype *> sampling_genotypes;
    std::vector<double> relative_infectivity_each_pp;

    for (int if_index = 0; if_index < interrupted_feeding_indices.size(); ++if_index) {
      // clear() is used to avoid memory reallocation
      sampling_genotypes.clear();
      relative_infectivity_each_pp.clear();

      if (config->within_host_induced_free_recombination()) {
        // get all infectious parasites from first person
        get_genotypes_profile_from_person(first_sampling[if_index], sampling_genotypes, relative_infectivity_each_pp);

        if (sampling_genotypes.empty()) {
          LOG(FATAL) << "first person has no infectious parasites, log10_total_infectious_denstiy = "
                     << first_sampling[if_index]->all_clonal_parasite_populations()->log10_total_infectious_denstiy;
        }

        if (interrupted_feeding_indices[if_index]) {
          auto temp_if = if_index;
          while (second_sampling[temp_if] == first_sampling[if_index]) {
            temp_if = random->random_uniform(second_sampling.size());
          }
          // interrupted feeding occurs
          get_genotypes_profile_from_person(second_sampling[temp_if], sampling_genotypes, relative_infectivity_each_pp);
          Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_count()[loc][res_23_list.size()]++;
        }

        if (sampling_genotypes.empty()) {
          LOG(FATAL) << "sampling_genotypes should not be empty";
        }
      } else {
        sampling_genotypes.clear();
        relative_infectivity_each_pp.clear();
        get_genotypes_profile_from_person(first_sampling[if_index], sampling_genotypes, relative_infectivity_each_pp);
        // get exactly 1 infectious parasite from first person
        auto first_genotype =
            random->roulette_sampling_tuple<Genotype>(1, relative_infectivity_each_pp, sampling_genotypes, false)[0];

        std::tuple<Genotype *, double> second_genotype = std::make_tuple(nullptr, 0.0);

        if (interrupted_feeding_indices[if_index]) {
          auto temp_if = if_index;
          while (second_sampling[temp_if] == first_sampling[if_index]) {
            temp_if = random->random_uniform(second_sampling.size());
          }
          sampling_genotypes.clear();
          relative_infectivity_each_pp.clear();
          get_genotypes_profile_from_person(second_sampling[temp_if], sampling_genotypes, relative_infectivity_each_pp);

          if (sampling_genotypes.size() > 0) {
            second_genotype = random->roulette_sampling_tuple<Genotype>(1, relative_infectivity_each_pp,
                                                                        sampling_genotypes, false)[0];
          }
          Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_count()[loc][res_23_list.size()]++;
        }

        sampling_genotypes.clear();
        relative_infectivity_each_pp.clear();
        sampling_genotypes.push_back(std::get<0>(first_genotype));
        relative_infectivity_each_pp.push_back(std::get<1>(first_genotype));

        if (std::get<0>(second_genotype) != nullptr) {
          sampling_genotypes.push_back(std::get<0>(second_genotype));
          relative_infectivity_each_pp.push_back(std::get<1>(second_genotype));
        }
      }

      auto parent_genotypes =
          random->roulette_sampling<Genotype>(2, relative_infectivity_each_pp, sampling_genotypes, false);

      Genotype *sampled_genotype =
          (parent_genotypes[0]->aa_sequence == parent_genotypes[1]->aa_sequence)
              ? parent_genotypes[0]
              : Genotype::free_recombine(config, random, parent_genotypes[0], parent_genotypes[1]);

      genotypes_table[tracking_index][loc][if_index] = sampled_genotype;

        for (int therapy_id = 6; therapy_id <= 8; therapy_id++) {
            auto *sc_therapy = dynamic_cast<SCTherapy *>(Model::CONFIG->therapy_db()[therapy_id]);
            if ((((parent_genotypes[0]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[0])) == drug_id_min_ec50[sc_therapy->drug_ids[0]])
                  && (parent_genotypes[0]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[1])) != drug_id_min_ec50[sc_therapy->drug_ids[1]])
                  && (parent_genotypes[1]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[0])) != drug_id_min_ec50[sc_therapy->drug_ids[0]])
                  && (parent_genotypes[1]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[1])) == drug_id_min_ec50[sc_therapy->drug_ids[1]]))
                 || ((parent_genotypes[1]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[0])) == drug_id_min_ec50[sc_therapy->drug_ids[0]])
                     && (parent_genotypes[1]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[1])) != drug_id_min_ec50[sc_therapy->drug_ids[1]])
                     && (parent_genotypes[0]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[0])) != drug_id_min_ec50[sc_therapy->drug_ids[0]])
                     && (parent_genotypes[0]->get_EC50_power_n(Model::CONFIG->drug_db()->at(sc_therapy->drug_ids[1])) == drug_id_min_ec50[sc_therapy->drug_ids[1]])))) {
                for(int res_id = 0; res_id < res_23_list.size(); res_id++){
                    if(get_resistant_strength_23(sampled_genotype,res_23_list[res_id])){
                        Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_count()[loc][res_id]++;
                    }
                }
                printf("\n");
            }
        }
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

bool Mosquito::string_contain(std::string str, std::string pattern) {
  if (str.find(pattern) != std::string::npos) {
    return true;
  }
  return false;
}

std::string Mosquito::get_resistant_strength_1(Genotype *genotype, std::string therapy) {
  std::string aa_seq = genotype->get_aa_sequence();
  std::vector<std::string> pattern_chromosome = split_string(aa_seq, '|');
  std::vector<std::string> chromosome_allele;
  int drug_num = 0;
  int allele_num = 0;
//  printf("Therapy: %s\n",therapy.c_str());
//  for(int i = 0; i < pattern_chromosome.size(); i++){
//    printf("%s|",pattern_chromosome[i].c_str());
//  }
//  printf("\n");
  for (int i = 0; i < pattern_chromosome.size(); i++) {
    std::string chromosome = "";
    for (int j = 0; j < pattern_chromosome[i].size(); j++) {
      chromosome += "0";
    }
    chromosome_allele.push_back(chromosome);
  }
  if (string_contain(split_string(therapy, '-')[0], "A")) //ART
  {
    if (pattern_chromosome[12].substr(10, 1) == "Y") {
      chromosome_allele[12].replace(10, 1, "1");
      drug_num += 1;
      allele_num += 1;
    }
  }
  if (string_contain(split_string(therapy, '-')[1], "L")) //LUM
  {
    if (pattern_chromosome[4].substr(0, 2) != "YY") {
      drug_num += 1;
    } else if (pattern_chromosome[6].substr(0, 1) != "T") {
      drug_num += 1;
    } else {
      drug_num = drug_num;
    }
    if (pattern_chromosome[6][0] == 'K' || pattern_chromosome[4][0] == 'N' || pattern_chromosome[4][1] == 'F' ||
        pattern_chromosome[4][2] == '2') {
      drug_num += 1;
    }
    if (drug_num > 0) {
      if (pattern_chromosome[6].substr(0, 1) == "K") { // K76T
        chromosome_allele[6].replace(0, 1, "1");
        allele_num += 1;
      }
      if (pattern_chromosome[4].substr(0, 1) == "N") {  // N86Y
        chromosome_allele[4].replace(0, 1, "1");
        allele_num += 1;
      }
      if (pattern_chromosome[4].substr(1, 1) == "F") {  // Y184F
        chromosome_allele[4].replace(1, 1, "1");
        allele_num += 1;
      }
    }
  }
  if (string_contain(split_string(therapy, '-')[1], "AQ")) //AQ
  {
    if (pattern_chromosome[6][0] == 'T' || pattern_chromosome[4][0] == 'Y' || pattern_chromosome[4][1] == 'Y') {
      drug_num += 1;
    }
    if (drug_num > 0) {
      if (pattern_chromosome[6].substr(0, 1) == "T") { // K76T
        chromosome_allele[6].replace(0, 1, "1");
        allele_num += 1;
      }
      if (pattern_chromosome[4].substr(0, 1) == "Y") {  // N86Y
        chromosome_allele[4].replace(0, 1, "1");
        allele_num += 1;
      }
      if (pattern_chromosome[4].substr(1, 1) == "Y") {  // Y184F
        chromosome_allele[4].replace(1, 1, "1");
        allele_num += 1;
      }
    }
  }
  if (string_contain(split_string(therapy, '-')[1], "PPQ")) //PPQ
  {
    if (pattern_chromosome[13].substr(0, 1) == "2") {
      chromosome_allele[13].replace(0, 1, "1");
      drug_num += 1;
      allele_num += 1;
    }
  }
  std::string pattern = std::to_string(drug_num) + "-" + std::to_string(allele_num);
//  printf("aa_seq: %s, pattern: %s, chromosome_allele:\n",aa_seq.c_str(),pattern.c_str());
//  for(int i = 0; i < pattern_chromosome.size(); i++){
//    printf("%s|",chromosome_allele[i].c_str());
//  }
  return pattern;
}

bool Mosquito::get_resistant_strength_23(Genotype *genotype, std::string resistance) {
  std::string aa_seq = genotype->get_aa_sequence();
  std::vector<std::string> pattern_chromosome = split_string(aa_seq, '|');
  std::vector<std::string> chromosome_allele;
  if (resistance == "AL:2") {
    if (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[6].substr(0, 1) == "K" \
         && pattern_chromosome[4].substr(0, 1) == "N" && pattern_chromosome[4].substr(1, 1) == "F" \
         && pattern_chromosome[13].substr(0, 1) == "1") {
      return true;
    }
  }

  if (resistance == "ASAQ:2") {
    if (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[6].substr(0, 1) == "T" \
         && pattern_chromosome[4].substr(0, 1) == "Y" && pattern_chromosome[4].substr(1, 1) == "Y" \
         && pattern_chromosome[13].substr(0, 1) == "1") {
      return true;
    }
  }

  if (resistance == "DHA-PPQ:2") {
    if (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[13].substr(0, 1) == "2") {
      return true;
    }
  }

  if (resistance == "DHA-PPQ-LUM:3") {
    if (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[6].substr(0, 1) == "K" \
         && pattern_chromosome[4].substr(0, 1) == "N" && pattern_chromosome[4].substr(1, 1) == "F" \
         && pattern_chromosome[13].substr(0, 1) == "2") {
      return true;
    }
  }

  if (resistance == "DHA-PPQ-AQ:3") {
    if (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[6].substr(0, 1) == "T" \
         && pattern_chromosome[4].substr(0, 1) == "Y" && pattern_chromosome[4].substr(1, 1) == "Y" \
         && pattern_chromosome[13].substr(0, 1) == "2") {
      return true;
    }
  }
  return false;
}
