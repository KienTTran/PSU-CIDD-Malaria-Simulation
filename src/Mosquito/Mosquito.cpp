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

      if (sampled_genotypes.size() == 1){
            Model::DATA_COLLECTOR->mosquito_sampled_genotypes_has_one_genotype()[loc]++;
            Model::DATA_COLLECTOR->monthly_mosquito_sampled_genotypes_has_one_genotype()[loc]++;
      }

      Genotype *sampled_genotype =
          (parent_genotypes[0]->aa_sequence == parent_genotypes[1]->aa_sequence)
              ? parent_genotypes[0]
              : Genotype::free_recombine(config, random, parent_genotypes[0], parent_genotypes[1]);

      genotypes_table[tracking_index][loc][if_index] = sampled_genotype;

      //Count DHA-PPQ(8) ASAQ(7) AL(6)
      //Count if male genotype resists to one drug and female genotype resists to another drug only, right now work on double resistant only
      //when genotype ec50_power_n == min_ec50, it is sensitive to that drug
        if (Model::SCHEDULER->current_time() >= Model::CONFIG->start_of_comparison_period())
        {
            for (int resistant_drug_pair_id = 0; resistant_drug_pair_id < resistant_drug_list.size(); resistant_drug_pair_id++) {
                auto drugs = resistant_drug_list[resistant_drug_pair_id].second;
                auto resistant_types = resistant_drug_list[resistant_drug_pair_id].first.size();
                for (int resistant_type_id = 0; resistant_type_id < resistant_types; resistant_type_id++) {
                    if(std::get<0>(count_no_condition_resistant_genotypes(parent_genotypes[0],
                                                             resistant_drug_pair_id,
                                                             resistant_type_id))){
                        Model::DATA_COLLECTOR->mosquito_sampled_male_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                        Model::DATA_COLLECTOR->monthly_mosquito_sampled_male_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                    }
                    if(std::get<0>(count_no_condition_resistant_genotypes(parent_genotypes[1],
                                                             resistant_drug_pair_id,
                                                             resistant_type_id))){
                        Model::DATA_COLLECTOR->mosquito_sampled_female_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                        Model::DATA_COLLECTOR->monthly_mosquito_sampled_female_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                    }
                    if (std::get<0>(count_one_condition_resistant_genotypes(config, loc, parent_genotypes, sampled_genotype, drugs,
                                              resistant_drug_pair_id,
                                              resistant_type_id,
                                              false))) {
                        Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                        Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                    }
                    if (std::get<0>(count_two_condition_resistant_genotypes(config, loc, parent_genotypes, sampled_genotype, drugs,
                                                                            resistant_drug_pair_id,
                                                                            resistant_type_id,
                                                                            true))) {
                        Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                        Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                    }
                    if(std::get<0>(count_no_condition_resistant_genotypes(sampled_genotype,
                                                             resistant_drug_pair_id,
                                                             resistant_type_id))){
                        Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                        Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]++;
                    }
                }
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
    if (Model::SCHEDULER->current_time() >= Model::CONFIG->start_of_comparison_period())
    {
      if(genotypes_table[tracking_index][location][genotype_index]->resistant_recombinations_in_mosquito.size() > 0){
          for(int res_drug_id = 0; res_drug_id < genotypes_table[tracking_index][location][genotype_index]->resistant_recombinations_in_mosquito.size(); res_drug_id++){
              MosquitoRecombinedGenotypeInfo resistant_info = genotypes_table[tracking_index][location][genotype_index]->resistant_recombinations_in_mosquito[res_drug_id];
              int res_drug_pair_id = resistant_info.second.first;
              int resistant_type_id = resistant_info.second.second;
              Model::DATA_COLLECTOR->mosquito_inflict_resistant_genotype_count()[location][res_drug_pair_id][resistant_type_id]++;
              Model::DATA_COLLECTOR->monthly_mosquito_inflict_resistant_genotype_count()[location][res_drug_pair_id][resistant_type_id]++;
          }
      }
    }
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

Mosquito::resistant_result_info Mosquito::count_no_condition_resistant_genotypes(Genotype *genotype, int resistant_drug_pair_id,int resistant_type_id) {
  std::vector<std::string> pattern_chromosome = split_string(genotype->aa_sequence, '|');
  std::vector<std::string> chromosome_allele;
    int res_points = 0;
    int mut_points = 0;
    std::string resistant_strength = "";
    bool is_resistant = false;
    //DHA-PPQ:2-2
    if (resistant_drug_pair_id == 0) {
        res_points = 2;
        mut_points = 2;
        if(resistant_type_id == 0) is_resistant = (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[13].substr(0, 1) == "2");//580Y-2
        resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
    }
    //ASAQ
    if (resistant_drug_pair_id == 1) {
        if(pattern_chromosome[12].substr(10, 1) == "Y"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[6].substr(0, 1) == "T"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(0, 1) == "Y"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(1, 1) == "Y"){
            mut_points++;
        }
        if(mut_points > 1){
            res_points++;
        }
        //Note that 2-4 resistant count includes 2-3 and 2-2 count
        if(resistant_type_id == 0) is_resistant = (res_points == 2 && mut_points == 2);//ASAQ:2-2
        if(resistant_type_id == 1) is_resistant = (res_points == 2 && mut_points == 3);//ASAQ:2-3
        if(resistant_type_id == 2) is_resistant = (res_points == 2 && mut_points == 4);//ASAQ:2-4
        resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
        if(resistant_type_id == 3) {
            is_resistant = (res_points == 2);//ASAQ:2
            resistant_strength = "";
        }
    }
    //AL
    if (resistant_drug_pair_id == 2) {
        if(pattern_chromosome[12].substr(10, 1) == "Y"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[6].substr(0, 1) == "K"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(0, 1) == "N"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(1, 1) == "F"){
            mut_points++;
        }
        if(mut_points > 1){
            res_points++;
        }
        if(resistant_type_id == 0) is_resistant = (res_points == 2 && mut_points == 2);//AL:2-2
        if(resistant_type_id == 1) is_resistant = (res_points == 2 && mut_points == 3);//AL:2-3
        if(resistant_type_id == 2) is_resistant = (res_points == 2 && mut_points == 4);//AL:2-4
        resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
        if(resistant_type_id == 3) {
            is_resistant = (res_points == 2);//AL:2
            resistant_strength = "";
        }
    }
    //DHA-PPQ-AQ
    if (resistant_drug_pair_id == 3) {
        if(pattern_chromosome[12].substr(10, 1) == "Y"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[13].substr(0, 1) == "2"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[6].substr(0, 1) == "T"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(0, 1) == "Y"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(1, 1) == "Y"){
            mut_points++;
        }
        if(mut_points > 2){
            res_points++;
        }
        if(resistant_type_id == 0) is_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-AQ:3-3
        if(resistant_type_id == 1) is_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-AQ:3-4
        if(resistant_type_id == 2) is_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-AQ:3-5
        resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
        if(resistant_type_id == 3) {
            is_resistant = (res_points == 3);//DHA-PPQ-AQ:3
            resistant_strength = "";
        }
    }
    //DHA-PPQ-LUM
    if (resistant_drug_pair_id == 4) {
        if(pattern_chromosome[12].substr(10, 1) == "Y"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[13].substr(0, 1) == "2"){
            res_points++;
            mut_points++;
        }
        if(pattern_chromosome[6].substr(0, 1) == "K"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(0, 1) == "N"){
            mut_points++;
        }
        if(pattern_chromosome[4].substr(1, 1) == "F"){
            mut_points++;
        }
        if(mut_points > 2){
            res_points++;
        }
        if(resistant_type_id == 0) is_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-LUM:3-3
        if(resistant_type_id == 1) is_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-LUM:3-4
        if(resistant_type_id == 2) is_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-LUM:3-5
        resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
        if(resistant_type_id == 3) {
            is_resistant = (res_points == 3);//DHA-PPQ-LUM:3
            resistant_strength = "";
        }
    }
    return std::make_tuple(is_resistant,resistant_drug_pair_id,resistant_type_id,resistant_strength);
}

Mosquito::resistant_result_info Mosquito::count_one_condition_resistant_genotypes(Config* config, int loc, std::vector<Genotype*> parent_genotypes, Genotype *genotype,
                                                                std::vector<int> drugs, int resistant_drug_pair_id, int resistant_type_id, bool verbose){
    //Double resistant
    if(drugs.size() == 2){
        if((parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))) //R0-R1
           ||(parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))) //R1-R0
           ){
            std::vector<std::string> pattern_chromosome = split_string(genotype->aa_sequence, '|');
            bool is_double_resistant = false;
            int res_points = 0;
            int mut_points = 0;
            std::string resistant_strength = "";
            //DHA-PPQ:2-2
            if (resistant_drug_pair_id == 0) {
                res_points = 2;
                mut_points = 2;
                if(resistant_type_id == 0) is_double_resistant = (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[13].substr(0, 1) == "2");//580Y-2
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
            }
            //ASAQ
            if (resistant_drug_pair_id == 1) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "T"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "Y"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "Y"){
                    mut_points++;
                }
                if(mut_points > 1){
                    res_points++;
                }
                //Note that 2-4 resistant count includes 2-3 and 2-2 count
                if(resistant_type_id == 0) is_double_resistant = (res_points == 2 && mut_points == 2);//ASAQ:2-2
                if(resistant_type_id == 1) is_double_resistant = (res_points == 2 && mut_points == 3);//ASAQ:2-3
                if(resistant_type_id == 2) is_double_resistant = (res_points == 2 && mut_points == 4);//ASAQ:2-4
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_double_resistant = (res_points == 2);//ASAQ:2
                    resistant_strength = "";
                }
            }
            //AL
            if (resistant_drug_pair_id == 2) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "K"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "N"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "F"){
                    mut_points++;
                }
                if(mut_points > 1){
                    res_points++;
                }
                if(resistant_type_id == 0) is_double_resistant = (res_points == 2 && mut_points == 2);//AL:2-2
                if(resistant_type_id == 1) is_double_resistant = (res_points == 2 && mut_points == 3);//AL:2-3
                if(resistant_type_id == 2) is_double_resistant = (res_points == 2 && mut_points == 4);//AL:2-4
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_double_resistant = (res_points == 2);//AL:2
                    resistant_strength = "";
                }
            }
            if(verbose && is_double_resistant){
                VLOG(1) << fmt::format("Count one condition {} resistant_drug_pair_id: {}\n"
                                       "genotype_m = \"{}\";\n"
                                       "genotype_f = \"{}\";\n"
                                       "genotype_c = \"{}\";\n"
                                       "m_ec50-d0: {:.10f}\tm_ec50-d1: {:.10f}\n"
                                       "f_ec50-d0: {:.10f}\tf_ec50-d1: {:.10f}\n"
                                       "min_ec50-d0: {:.10f}\tmin_ec50-d1: {:.10f}\n"
                                       "resistant_type: {} {} {}",
                                       Model::SCHEDULER->current_time(),
                                       resistant_drug_pair_id,
                                       parent_genotypes[0]->get_aa_sequence().c_str(),
                                       parent_genotypes[1]->get_aa_sequence().c_str(),
                                       genotype->aa_sequence.c_str(),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       drug_id_min_ec50[drugs[0]],
                                       drug_id_min_ec50[drugs[1]],
                                       resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                       resistant_strength,
                                       is_double_resistant);
                VLOG(1) << fmt::format("Count one condition [{}][{}][{}]: month_clonal_resistant: {}\tmonth_mos_resistant: {}\tmonth_mos_inflict: {}\tcumm_clonal_resistant: {}\tcumm_mos_resistant: {}\tcumm_mos_inflict: {}",
                                       loc,resistant_drug_pair_id,resistant_type_id,
                                       Model::DATA_COLLECTOR->monthly_clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->monthly_mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]);
            }
            return std::make_tuple(is_double_resistant,resistant_drug_pair_id,resistant_type_id,resistant_strength);
        }
    }
    //Triple resistant
    if(drugs.size() == 3){
        if(((parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))    && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))) //0-12
        || (parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))) //1-02
        || (parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))     && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))))//2-01
        || ((parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))      && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))) //0-12 g0 <-> g1
        || (parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))) //0-12
        || (parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))     && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))))//1-02
        ){
            std::vector<std::string> pattern_chromosome = split_string(genotype->aa_sequence, '|');
            bool is_triple_resistant = false;
            int res_points = 0;
            int mut_points = 0;
            std::string resistant_strength = "";
            //DHA-PPQ-AQ
            if (resistant_drug_pair_id == 3) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[13].substr(0, 1) == "2"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "T"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "Y"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "Y"){
                    mut_points++;
                }
                if(mut_points > 2){
                    res_points++;
                }
                if(resistant_type_id == 0) is_triple_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-AQ:3-3
                if(resistant_type_id == 1) is_triple_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-AQ:3-4
                if(resistant_type_id == 2) is_triple_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-AQ:3-5
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_triple_resistant = (res_points == 3);//DHA-PPQ-AQ:3
                    resistant_strength = "";
                }
            }
            //DHA-PPQ-LUM
            if (resistant_drug_pair_id == 4) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[13].substr(0, 1) == "2"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "K"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "N"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "F"){
                    mut_points++;
                }
                if(mut_points > 2){
                    res_points++;
                }
                if(resistant_type_id == 0) is_triple_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-LUM:3-3
                if(resistant_type_id == 1) is_triple_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-LUM:3-4
                if(resistant_type_id == 2) is_triple_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-LUM:3-5
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_triple_resistant = (res_points == 3);//DHA-PPQ-LUM:3
                    resistant_strength = "";
                }
            }
            if(verbose && is_triple_resistant){
                VLOG(1) << fmt::format("Count one condition {} resistant_drug_pair_id: {} \n"
                                       "genotype_m = \"{}\";\n"
                                       "genotype_f = \"{}\";\n"
                                       "genotype_c = \"{}\";\n"
                                       "m_ec50-d0: {:.10f}\tm_ec50-d1: {:.10f}\tm_ec50-d2: {:.10f}\n"
                                       "f_ec50-d0: {:.10f}\tf_ec50-d1: {:.10f}\tf_ec50-d2: {:.10f}\n"
                                       "min_ec50-d0: {:.10f}\tmin_ec50-d1: {:.10f}\tmin_ec50-d2: {:.10f}\n"
                                       "resistant_type: {} {} {}",
                                       Model::SCHEDULER->current_time(),
                                       resistant_drug_pair_id,
                                       parent_genotypes[0]->get_aa_sequence().c_str(),
                                       parent_genotypes[1]->get_aa_sequence().c_str(),
                                       genotype->aa_sequence.c_str(),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[2])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[2])),
                                       drug_id_min_ec50[drugs[0]],
                                       drug_id_min_ec50[drugs[1]],
                                       drug_id_min_ec50[drugs[2]],
                                       resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                       resistant_strength,
                                       is_triple_resistant);
                VLOG(1) << fmt::format("Count one condition [{}][{}][{}]: month_clonal_resistant: {}\tmonth_mos_resistant: {}\tmonth_mos_inflict: {}\tcumm_clonal_resistant: {}\tcumm_mos_resistant: {}\tcumm_mos_inflict: {}",
                                       loc,resistant_drug_pair_id,resistant_type_id,
                                       Model::DATA_COLLECTOR->monthly_clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->monthly_mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_one_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]);
            }
            return std::make_tuple(is_triple_resistant,resistant_drug_pair_id,resistant_type_id,resistant_strength);
        }
    }
  return std::make_tuple(false,resistant_drug_pair_id,resistant_type_id,"0-0");
}

Mosquito::resistant_result_info Mosquito::count_two_condition_resistant_genotypes(Config* config, int loc, std::vector<Genotype*> parent_genotypes, Genotype *genotype,
                                                                    std::vector<int> drugs, int resistant_drug_pair_id, int resistant_type_id, bool verbose){
    //Double resistant
    if(drugs.size() == 2){
        if((parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1])) //g0 - g1
        && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))   && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))) //RS-SR
        ||(parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1])) //g1 - g0
        && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))    && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))) //RS-SR
        ){
            std::vector<std::string> pattern_chromosome = split_string(genotype->aa_sequence, '|');
            bool is_double_resistant = false;
            int res_points = 0;
            int mut_points = 0;
            std::string resistant_strength = "";
            //DHA-PPQ:2-2
            if (resistant_drug_pair_id == 0) {
                res_points = 2;
                mut_points = 2;
                if(resistant_type_id == 0) is_double_resistant = (pattern_chromosome[12].substr(10, 1) == "Y" && pattern_chromosome[13].substr(0, 1) == "2");//580Y-2
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
            }
            //ASAQ
            if (resistant_drug_pair_id == 1) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "T"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "Y"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "Y"){
                    mut_points++;
                }
                if(mut_points > 1){
                    res_points++;
                }
                //Note that 2-4 resistant count includes 2-3 and 2-2 count
                if(resistant_type_id == 0) is_double_resistant = (res_points == 2 && mut_points == 2);//ASAQ:2-2
                if(resistant_type_id == 1) is_double_resistant = (res_points == 2 && mut_points == 3);//ASAQ:2-3
                if(resistant_type_id == 2) is_double_resistant = (res_points == 2 && mut_points == 4);//ASAQ:2-4
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_double_resistant = (res_points == 2);//ASAQ:2
                    resistant_strength = "";
                }
            }
            //AL
            if (resistant_drug_pair_id == 2) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "K"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "N"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "F"){
                    mut_points++;
                }
                if(mut_points > 1){
                    res_points++;
                }
                if(resistant_type_id == 0) is_double_resistant = (res_points == 2 && mut_points == 2);//AL:2-2
                if(resistant_type_id == 1) is_double_resistant = (res_points == 2 && mut_points == 3);//AL:2-3
                if(resistant_type_id == 2) is_double_resistant = (res_points == 2 && mut_points == 4);//AL:2-4
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_double_resistant = (res_points == 2);//AL:2
                    resistant_strength = "";
                }
            }
            if(is_double_resistant){
                std::vector<std::pair<int,std::string>> parent_info;
                parent_info.push_back(std::make_pair(parent_genotypes[0]->genotype_id,parent_genotypes[0]->aa_sequence));
                parent_info.push_back(std::make_pair(parent_genotypes[1]->genotype_id,parent_genotypes[1]->aa_sequence));
                MosquitoRecombinedGenotypeInfo resistant_info = std::make_pair(parent_info,std::make_pair(resistant_drug_pair_id,resistant_type_id));
                genotype->resistant_recombinations_in_mosquito.push_back(resistant_info);
                parent_info.clear();
            }
            if(verbose && is_double_resistant){
                VLOG(1) << fmt::format("Count two condition {} resistant_drug_pair_id: {}\n"
                                       "genotype_m = \"{}\";\n"
                                       "genotype_f = \"{}\";\n"
                                       "genotype_c = \"{}\";\n"
                                       "m_ec50-d0: {:.10f}\tm_ec50-d1: {:.10f}\n"
                                       "f_ec50-d0: {:.10f}\tf_ec50-d1: {:.10f}\n"
                                       "min_ec50-d0: {:.10f}\tmin_ec50-d1: {:.10f}\n"
                                       "resistant_type: {} {} {}",
                                       Model::SCHEDULER->current_time(),
                                       resistant_drug_pair_id,
                                       parent_genotypes[0]->get_aa_sequence().c_str(),
                                       parent_genotypes[1]->get_aa_sequence().c_str(),
                                       genotype->aa_sequence.c_str(),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       drug_id_min_ec50[drugs[0]],
                                       drug_id_min_ec50[drugs[1]],
                                       resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                       resistant_strength,
                                       is_double_resistant);
                VLOG(1) << fmt::format("Count two condition [{}][{}][{}]: month_clonal_resistant: {}\tmonth_mos_resistant: {}\tmonth_mos_inflict: {}\tcumm_clonal_resistant: {}\tcumm_mos_resistant: {}\tcumm_mos_inflict: {}",
                                       loc,resistant_drug_pair_id,resistant_type_id,
                                       Model::DATA_COLLECTOR->monthly_clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->monthly_mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]);
            }
            return std::make_tuple(is_double_resistant,resistant_drug_pair_id,resistant_type_id,resistant_strength);
        }
    }
    //Triple resistant
    if(drugs.size() == 3){
        if(((parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))    && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))     && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))//g0 - g1
        && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))     && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))) //RSS-SRR
        || (!parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))    && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))      && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))
        && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))      && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))) //SRS-RSR
        || (!parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))    && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))
        && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))      && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2])))) //SSR-RRS
        || ((parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))     && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2])) //g1 - g0
        && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))     && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))      && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))) //RSS-SRR
        || (!parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))      && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))
        && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2]))) //SRS-RSR
        || (!parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[0]))    && !parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[1]))     && parent_genotypes[1]->resist_to(config->drug_db()->at(drugs[2]))
        && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[0]))      && parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[1]))      && !parent_genotypes[0]->resist_to(config->drug_db()->at(drugs[2])))) //SSR-RRS
        ){
            std::vector<std::string> pattern_chromosome = split_string(genotype->aa_sequence, '|');
            bool is_triple_resistant = false;
            int res_points = 0;
            int mut_points = 0;
            std::string resistant_strength = "";
            //DHA-PPQ-AQ
            if (resistant_drug_pair_id == 3) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[13].substr(0, 1) == "2"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "T"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "Y"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "Y"){
                    mut_points++;
                }
                if(mut_points > 2){
                    res_points++;
                }
                if(resistant_type_id == 0) is_triple_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-AQ:3-3
                if(resistant_type_id == 1) is_triple_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-AQ:3-4
                if(resistant_type_id == 2) is_triple_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-AQ:3-5
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_triple_resistant = (res_points == 3);//DHA-PPQ-AQ:3
                    resistant_strength = "";
                }
            }
            //DHA-PPQ-LUM
            if (resistant_drug_pair_id == 4) {
                if(pattern_chromosome[12].substr(10, 1) == "Y"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[13].substr(0, 1) == "2"){
                    res_points++;
                    mut_points++;
                }
                if(pattern_chromosome[6].substr(0, 1) == "K"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(0, 1) == "N"){
                    mut_points++;
                }
                if(pattern_chromosome[4].substr(1, 1) == "F"){
                    mut_points++;
                }
                if(mut_points > 2){
                    res_points++;
                }
                if(resistant_type_id == 0) is_triple_resistant = (res_points == 3 && mut_points == 3);//DHA-PPQ-LUM:3-3
                if(resistant_type_id == 1) is_triple_resistant = (res_points == 3 && mut_points == 4);//DHA-PPQ-LUM:3-4
                if(resistant_type_id == 2) is_triple_resistant = (res_points == 3 && mut_points == 5);//DHA-PPQ-LUM:3-5
                resistant_strength = std::to_string(res_points) + "-" + std::to_string(mut_points);
                if(resistant_type_id == 3) {
                    is_triple_resistant = (res_points == 3);//DHA-PPQ-LUM:3
                    resistant_strength = "";
                }
            }
            if(is_triple_resistant){
                std::vector<std::pair<int,std::string>> parent_info;
                parent_info.push_back(std::make_pair(parent_genotypes[0]->genotype_id,parent_genotypes[0]->aa_sequence));
                parent_info.push_back(std::make_pair(parent_genotypes[1]->genotype_id,parent_genotypes[1]->aa_sequence));
                MosquitoRecombinedGenotypeInfo resistant_info = std::make_pair(parent_info,std::make_pair(resistant_drug_pair_id,resistant_type_id));
                genotype->resistant_recombinations_in_mosquito.push_back(resistant_info);
                parent_info.clear();
            }
            if(verbose && is_triple_resistant){
                VLOG(1) << fmt::format("Count two condition {} resistant_drug_pair_id: {} \n"
                                       "genotype_m = \"{}\";\n"
                                       "genotype_f = \"{}\";\n"
                                       "genotype_c = \"{}\";\n"
                                       "m_ec50-d0: {:.10f}\tm_ec50-d1: {:.10f}\tm_ec50-d2: {:.10f}\n"
                                       "f_ec50-d0: {:.10f}\tf_ec50-d1: {:.10f}\tf_ec50-d2: {:.10f}\n"
                                       "min_ec50-d0: {:.10f}\tmin_ec50-d1: {:.10f}\tmin_ec50-d2: {:.10f}\n"
                                       "resistant_type: {} {} {}",
                                       Model::SCHEDULER->current_time(),
                                       resistant_drug_pair_id,
                                       parent_genotypes[0]->get_aa_sequence().c_str(),
                                       parent_genotypes[1]->get_aa_sequence().c_str(),
                                       genotype->aa_sequence.c_str(),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[0]->get_EC50_power_n(config->drug_db()->at(drugs[2])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[0])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[1])),
                                       parent_genotypes[1]->get_EC50_power_n(config->drug_db()->at(drugs[2])),
                                       drug_id_min_ec50[drugs[0]],
                                       drug_id_min_ec50[drugs[1]],
                                       drug_id_min_ec50[drugs[2]],
                                       resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                       resistant_strength,
                                       is_triple_resistant);
                VLOG(1) << fmt::format("Count two condition [{}][{}][{}]: month_clonal_resistant: {}\tmonth_mos_resistant: {}\tmonth_mos_inflict: {}\tcumm_clonal_resistant: {}\tcumm_mos_resistant: {}\tcumm_mos_inflict: {}",
                                       loc,resistant_drug_pair_id,resistant_type_id,
                                       Model::DATA_COLLECTOR->monthly_clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->monthly_mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->monthly_mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->clonal_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id],
                                       Model::DATA_COLLECTOR->mosquito_recombined_resistant_genotype_two_condition_count()[loc][resistant_drug_pair_id][resistant_type_id] + 1,
                                       Model::DATA_COLLECTOR->mosquito_inflict_resistant_genotype_count()[loc][resistant_drug_pair_id][resistant_type_id]);
            }
            return std::make_tuple(is_triple_resistant,resistant_drug_pair_id,resistant_type_id,resistant_strength);
        }
    }
    return std::make_tuple(false,resistant_drug_pair_id,resistant_type_id,"0-0");
}

