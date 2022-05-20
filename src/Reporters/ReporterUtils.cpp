//
// Created by Nguyen Tran on 2019-04-18.
//

#include "ReporterUtils.h"

#include <Core/Config/Config.h>

#include <vector>

#include "Model.h"
#include "Population/ClonalParasitePopulation.h"
#include "Population/Properties/PersonIndexByLocationStateAgeClass.h"
#include "Population/SingleHostClonalParasitePopulations.h"

const std::string group_sep = "-1111\t";
const std::string sep = "\t";

void ReporterUtils::output_genotype_frequency1(std::stringstream& ss, const int& number_of_genotypes,
                                               PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  std::vector<double> result1_all(number_of_genotypes, 0.0);

  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result1(number_of_genotypes, 0.0);
    auto sum1 = 0.0;

    for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
      for (auto ac = 0; ac < number_of_age_classes; ac++) {
        const auto size = pi->vPerson()[loc][hs][ac].size();
        for (auto i = 0ull; i < size; i++) {
          auto* person = pi->vPerson()[loc][hs][ac][i];

          if (!person->all_clonal_parasite_populations()->parasites()->empty()) {
            sum1 += 1;
            sum1_all += 1;
          }
          // sum2 += person->all_clonal_parasite_populations()->parasites()->size();
          std::map<int, int> individual_genotype_map;

          for (auto* parasite_population : *(person->all_clonal_parasite_populations()->parasites())) {
            const auto g_id = parasite_population->genotype()->genotype_id;
            // result2[g_id] += 1;
            if (individual_genotype_map.find(g_id) == individual_genotype_map.end()) {
              individual_genotype_map[parasite_population->genotype()->genotype_id] = 1;
            } else {
              individual_genotype_map[parasite_population->genotype()->genotype_id] += 1;
            }
          }

          for (const auto genotype : individual_genotype_map) {
            result1[genotype.first] += 1;
            result1_all[genotype.first] += 1;
          }
        }
      }
    }

    for (auto& i : result1) {
      i /= sum1;
      ss << (sum1 == 0 ? 0 : i) << sep;
    }
  }
  ss << group_sep;
  for (auto& i : result1_all) {
    i /= sum1_all;
    ss << (sum1_all == 0 ? 0 : i) << sep;
  }
}

void ReporterUtils::output_genotype_frequency2(std::stringstream& ss, const int& number_of_genotypes,
                                               PersonIndexByLocationStateAgeClass* pi) {
  auto sum2_all = 0.0;
  std::vector<double> result2_all(number_of_genotypes, 0.0);
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result2(number_of_genotypes, 0.0);
    auto sum2 = 0.0;

    for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
      for (auto ac = 0; ac < number_of_age_classes; ac++) {
        const auto size = pi->vPerson()[loc][hs][ac].size();
        for (auto i = 0ull; i < size; i++) {
          auto* person = pi->vPerson()[loc][hs][ac][i];

          sum2 += person->all_clonal_parasite_populations()->parasites()->size();
          sum2_all += sum2;
          std::map<int, int> individual_genotype_map;

          for (auto* parasite_population : *(person->all_clonal_parasite_populations()->parasites())) {
            const auto g_id = parasite_population->genotype()->genotype_id;
            result2[g_id] += 1;
            result2_all[g_id] += 1;
          }
        }
      }
    }

    // output for each location
    for (auto& i : result2) {
      i /= sum2;
      ss << (sum2 == 0 ? 0 : i) << sep;
    }
  }
  ss << group_sep;
  // output for all locations
  for (auto& i : result2_all) {
    i /= sum2_all;
    ss << (sum2_all == 0 ? 0 : i) << sep;
  }
}

void ReporterUtils::output_genotype_frequency3(std::stringstream& ss, const int& number_of_genotypes,
                                               PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  std::vector<double> result3_all(number_of_genotypes, 0.0);
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result3(number_of_genotypes, 0.0);
    auto sum1 = 0.0;

    for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
      for (auto ac = 0; ac < number_of_age_classes; ac++) {
        const auto size = pi->vPerson()[loc][hs][ac].size();
        for (auto i = 0ull; i < size; i++) {
          auto* person = pi->vPerson()[loc][hs][ac][i];

          if (!person->all_clonal_parasite_populations()->parasites()->empty()) {
            sum1 += 1;
            sum1_all += 1;
          }

          std::map<int, int> individual_genotype_map;

          for (auto* parasite_population : *(person->all_clonal_parasite_populations()->parasites())) {
            const auto g_id = parasite_population->genotype()->genotype_id;
            if (individual_genotype_map.find(g_id) == individual_genotype_map.end()) {
              individual_genotype_map[parasite_population->genotype()->genotype_id] = 1;
            } else {
              individual_genotype_map[parasite_population->genotype()->genotype_id] += 1;
            }
          }

          for (const auto genotype : individual_genotype_map) {
            result3[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
            result3_all[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
          }
        }
      }
    }
    // output per location
    // TODO: implement dynamic way to output for each location
//
//    for (auto& i : result3) {
//      i /= sum1;
//      ss << (sum1 == 0 ? 0 : i) << sep;
//    }
  }
//  ss << group_sep;

  // this is for all locations
  for (auto& i : result3_all) {
    i /= sum1_all;
    ss << (sum1_all == 0 ? 0 : i) << sep;
  }
//
//  ss << group_sep;
//  ss << sum1_all << sep;
}

void ReporterUtils::output_3_genotype_frequency(std::stringstream& ss, const int& number_of_genotypes,
                                                PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  auto sum2_all = 0.0;
  std::vector<double> result1_all(number_of_genotypes, 0.0);
  std::vector<double> result2_all(number_of_genotypes, 0.0);
  std::vector<double> result3_all(number_of_genotypes, 0.0);

  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result1(number_of_genotypes, 0.0);
    std::vector<double> result2(number_of_genotypes, 0.0);
    std::vector<double> result3(number_of_genotypes, 0.0);

    auto sum2 = 0.0;
    auto sum1 = 0.0;

    for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
      for (auto ac = 0; ac < number_of_age_classes; ac++) {
        const auto size = pi->vPerson()[loc][hs][ac].size();
        for (auto i = 0ull; i < size; i++) {
          auto* person = pi->vPerson()[loc][hs][ac][i];

          if (!person->all_clonal_parasite_populations()->parasites()->empty()) {
            sum1 += 1;
            sum1_all += 1;
          }
          sum2 += person->all_clonal_parasite_populations()->parasites()->size();
          sum2_all += sum2;
          std::map<int, int> individual_genotype_map;

          for (auto* parasite_population : *(person->all_clonal_parasite_populations()->parasites())) {
            std::cout << "hello" << std::endl;
            const auto g_id = parasite_population->genotype()->genotype_id;
            result2[g_id] += 1;
            result2_all[g_id] += 1;
            if (individual_genotype_map.find(g_id) == individual_genotype_map.end()) {
              individual_genotype_map[parasite_population->genotype()->genotype_id] = 1;
            } else {
              individual_genotype_map[parasite_population->genotype()->genotype_id] += 1;
            }
          }

          for (const auto genotype : individual_genotype_map) {
            result1[genotype.first] += 1;
            result1_all[genotype.first] += 1;
            result3[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
            result3_all[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
          }
        }
      }
    }

    for (auto& i : result1) {
      i /= sum1;
    }

    for (auto& i : result2) {
      i /= sum2;
    }

    for (auto& i : result3) {
      i /= sum1;
    }

    for (auto j = 0; j < number_of_genotypes; ++j) {
      if (sum1 == 0) {
        ss << 0 << sep << 0 << sep << 0 << sep;
      } else {
        ss << result1[j] << sep;
        ss << result2[j] << sep;
        ss << result3[j] << sep;
      }
    }

    ss << group_sep;
  }

  // this is for all locations
  ss << group_sep;
  for (auto j = 0; j < number_of_genotypes; ++j) {
    if (sum1_all == 0) {
      ss << 0 << sep << 0 << sep << 0 << sep;
    } else {
      result1_all[j] /= sum1_all;
      result2_all[j] /= sum2_all;
      result3_all[j] /= sum1_all;

      ss << result1_all[j] << sep;
      ss << result1_all[j] << sep;
      ss << result1_all[j] << sep;
    }
  }
}

void ReporterUtils::output_moi(std::stringstream& ss, PersonIndexByLocationStateAgeClass* pi) {
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0ul; loc < number_of_locations; loc++) {
    auto pop_sum_location = 0;
    for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
      for (auto ac = 0ul; ac < number_of_age_classes; ac++) {
        std::size_t size = pi->vPerson()[loc][hs][ac].size();
        for (int i = 0; i < size; i++) {
          Person* p = pi->vPerson()[loc][hs][ac][i];
          int moi = p->all_clonal_parasite_populations()->size();
          ss << moi << "\n";
        }
      }
    }
  }
  CLOG(INFO, "moi_reporter") << ss.str();
  ss.str("");
}

void ReporterUtils::initialize_moi_file_logger() {
  const std::string OUTPUT_FORMAT = "[%level] [%logger] [%host] [%func] [%loc] %msg";

  el::Configurations moi_reporter_logger;
  moi_reporter_logger.setToDefault();
  moi_reporter_logger.set(el::Level::Debug, el::ConfigurationType::Format, OUTPUT_FORMAT);
  moi_reporter_logger.set(el::Level::Error, el::ConfigurationType::Format, OUTPUT_FORMAT);
  moi_reporter_logger.set(el::Level::Fatal, el::ConfigurationType::Format, OUTPUT_FORMAT);
  moi_reporter_logger.set(el::Level::Trace, el::ConfigurationType::Format, OUTPUT_FORMAT);
  moi_reporter_logger.set(el::Level::Info, el::ConfigurationType::Format, "%msg");
  moi_reporter_logger.set(el::Level::Warning, el::ConfigurationType::Format, "[%level] [%logger] %msg");
  moi_reporter_logger.set(el::Level::Verbose, el::ConfigurationType::Format, "[%level-%vlevel] [%logger] %msg");

  moi_reporter_logger.setGlobally(el::ConfigurationType::ToFile, "true");
  moi_reporter_logger.setGlobally(el::ConfigurationType::Filename,
                                  fmt::format("{}moi_{}.txt", "", Model::MODEL->cluster_job_number()));
  moi_reporter_logger.setGlobally(el::ConfigurationType::ToStandardOutput, "false");
  moi_reporter_logger.setGlobally(el::ConfigurationType::LogFlushThreshold, "100");
  // default logger uses default configurations
  el::Loggers::reconfigureLogger("moi_reporter", moi_reporter_logger);
}

//
//
//
// void MonthlyReporter::print_monthly_incidence_by_location() {
//  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
//    ss << Model::DATA_COLLECTOR->monthly_number_of_treatment_by_location()[loc] << sep;
//  }
//
//  ss << group_sep;
//
//  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
//    ss << Model::DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] << sep;
//  }
//}

inline std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

static std::map<std::string,int> mutated_allele_pos;
std::map<std::string,double> alleles_freq = std::map<std::string,double>();
std::map<std::string,double> genotype_alleles_freq = std::map<std::string,double>();

void ReporterUtils::output_genotype_frequency4(std::stringstream& ss, const int& number_of_genotypes,
                                               PersonIndexByLocationStateAgeClass* pi) {
    auto sum1_all = 0.0;
    std::vector<double> result3_all(number_of_genotypes, 0.0);
    const auto number_of_locations = pi->vPerson().size();
    const auto number_of_age_classes = pi->vPerson()[0][0].size();

    std::map<std::string,int> individual_with_mutated_alleles_count;

    for (auto loc = 0; loc < number_of_locations; loc++) {
        std::vector<double> result3(number_of_genotypes, 0.0);
        auto sum1 = 0.0;

        for (auto hs = 0; hs < Person::NUMBER_OF_STATE - 1; hs++) {
            for (auto ac = 0; ac < number_of_age_classes; ac++) {
                const auto size = pi->vPerson()[loc][hs][ac].size();
                for (auto i = 0ull; i < size; i++) {
                    auto* person = pi->vPerson()[loc][hs][ac][i];

                    if (!person->all_clonal_parasite_populations()->parasites()->empty()) {
                        sum1 += 1;
                        sum1_all += 1;
                    }

                    std::map<int, int> individual_genotype_map;

                    for (auto* parasite_population : *(person->all_clonal_parasite_populations()->parasites())) {

                        const auto person_genotype = parasite_population->genotype();
                        const auto g_id = person_genotype->genotype_id;

                        if (individual_genotype_map.find(g_id) == individual_genotype_map.end()) {
                            individual_genotype_map[parasite_population->genotype()->genotype_id] = 1;
                        } else {
                            individual_genotype_map[parasite_population->genotype()->genotype_id] += 1;
                        }
                    }

                    for (const auto genotype_stat : individual_genotype_map) {
                        result3[genotype_stat.first] +=
                                genotype_stat.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
                        result3_all[genotype_stat.first] +=
                                genotype_stat.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());

//                        const auto run_genotype = Model::CONFIG->genotype_db.at(genotype_stat.first);
//                        const auto g_id = run_genotype->genotype_id;
//
//                        if(g_id >= Model::CONFIG->initial_parasite_info().size()){
//                            for(int j = 0; j < Model::CONFIG->initial_parasite_info().size(); j++){
//                                auto init_genotype = Model::CONFIG->genotype_db.at(Model::CONFIG->initial_parasite_info()[j].parasite_type_id);
//                                for(int k = 0; k < run_genotype->get_aa_sequence().length(); k++) {
//                                    if(run_genotype->get_aa_sequence().at(k) != init_genotype->get_aa_sequence().at(k)){
//                                        std::string key = std::to_string(g_id);
//                                        key.push_back('-');
//                                        key.push_back(init_genotype->get_aa_sequence().at(k));
//                                        key.push_back('-');
//                                        key.push_back(run_genotype->get_aa_sequence().at(k));
//                                        if (individual_with_mutated_alleles_count.find(key) == individual_with_mutated_alleles_count.end()) {
//                                            individual_with_mutated_alleles_count[key] = 1;
//                                        } else {
//                                            individual_with_mutated_alleles_count[key] += 1;
//                                        }
//                                        mutated_allele_pos[key] = k;
//                                    }
//                                }
//                            }
//                        }
                    }
                }
            }
        }
        // output per location
        // TODO: implement dynamic way to output for each location
//
//    for (auto& i : result3) {
//      i /= sum1;
//      ss << (sum1 == 0 ? 0 : i) << sep;
//    }
    }
//  ss << group_sep;

    // this is for all locations
//    for (auto& i : result3_all) {
//        i /= sum1_all;
//        ss << (sum1_all == 0 ? 0 : i) << sep;
//    }

//    for(const auto allele_pos : mutated_allele_pos){
//        printf("Day %d Alelle %s pos = %d\n",Model::SCHEDULER->current_time(), allele_pos.first.c_str(), allele_pos.second);
//    }
//    for(const auto allele_count : individual_with_mutated_alleles_count){
//        printf("Day %d Alelle %s count = %d\n",Model::SCHEDULER->current_time(), allele_count.first.c_str(), allele_count.second);
//    }

    alleles_freq = std::map<std::string,double>();
    genotype_alleles_freq = std::map<std::string,double>();
    std::map<std::string,std::string> alleles;
    for(int i = 0; i < result3_all.size(); i++){
        for(int j = 0; j < Model::CONFIG->mutation_mask().length(); j++){
            if(Model::CONFIG->mutation_mask()[j] == '1'){
                alleles[std::to_string(i) + "-" + std::to_string(j)] = Model::CONFIG->genotype_db.at(i)->get_aa_sequence().substr(j,1);
            }
        }
    }

    for(const auto allele : alleles){
        printf("Day %d Allele loci %s is %s\n",Model::SCHEDULER->current_time(), allele.first.c_str(), allele.second.c_str());
    }

    for(int i = 0; i < result3_all.size(); i++){
        auto genotype = Model::CONFIG->genotype_db.at(i);
        double genotype_freq = 0.0;
        if(sum1_all == 0){
            genotype_freq = 0.0;
        }
        else{
            genotype_freq = result3_all[i] / sum1_all;
        }
        printf("Day %d Genotype id %d %s freq %f\n",Model::SCHEDULER->current_time(),i,Model::CONFIG->genotype_db.at(i)->get_aa_sequence().c_str(),genotype_freq);
//        for(auto map : mutated_allele_pos){
//            if(i == map.first.front() - 48){
//                if (alleles_freq.find(map.first) == alleles_freq.end()) {
//                    alleles_freq[map.first] = genotype_freq;
//                } else {
//                    alleles_freq[map.first] += genotype_freq;
//                }
//            }
//        }
        std::string genotype_aa = Model::CONFIG->genotype_db.at(i)->get_aa_sequence();
        for(const auto allele : alleles){
            std::string key = allele.first + "-" + allele.second;
            if(split(allele.first,'-').front() == std::to_string(i)){
                alleles_freq[key] = genotype_freq;
            }
        }
    }

    std::vector<std::string> added_genotype_id_mask_pos;
    std::vector<std::string> added_mask_pos;
    for(const auto freq : alleles_freq){
        std::string genotype_id = split(freq.first,'-').front();
        std::string mask_pos = split(freq.first, '-')[1];
        std::string allele = split(freq.first,'-').back();
        genotype_alleles_freq[mask_pos+"-"+allele] += freq.second;
    }

    for(const auto freq : alleles_freq){
        printf("Day %d Allele %s freq = %f\n",Model::SCHEDULER->current_time(), freq.first.c_str(), freq.second);
    }

    for(const auto freq : genotype_alleles_freq){
        printf("Day %d Genotype allele %s freq = %f\n",Model::SCHEDULER->current_time(), freq.first.c_str(), freq.second);
        ss << freq.second << sep;
    }

    if(alleles_freq.size() == 0) ss << 0 << sep;

    std::cout << ss.str() << std::endl;
//
//  ss << group_sep;
//  ss << sum1_all << sep;
}
