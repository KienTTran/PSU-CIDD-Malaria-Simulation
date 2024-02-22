//
// Created by Nguyen Tran on 2019-04-18.
//

#include "ReporterUtils.cuh"

#include <Core/Config/Config.h>

#include <vector>

#include "Model.h"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Mosquito/Mosquito.cuh"

const std::string group_sep = "-1111\t";
const std::string sep = "\t";

void GPU::ReporterUtils::output_genotype_frequency1(std::stringstream& ss, const int& number_of_genotypes,
                                               GPU::PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  std::vector<double> result1_all(number_of_genotypes, 0.0);

  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result1(number_of_genotypes, 0.0);
    auto sum1 = 0.0;

    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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

void GPU::ReporterUtils::output_genotype_frequency2(std::stringstream& ss, const int& number_of_genotypes,
                                               GPU::PersonIndexByLocationStateAgeClass* pi) {
  auto sum2_all = 0.0;
  std::vector<double> result2_all(number_of_genotypes, 0.0);
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result2(number_of_genotypes, 0.0);
    auto sum2 = 0.0;

    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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


void GPU::ReporterUtils::output_genotype_frequency3(std::stringstream& ss, const int& number_of_genotypes,
                                               GPU::PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  std::vector<double> result3_all(number_of_genotypes, 0.0);
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result3(number_of_genotypes, 0.0);
    auto sum1 = 0.0;

    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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
//            individual_genotype_cr[parasite_population->genotype()->genotype_id] = parasite_population->genotype()->daily_fitness_multiple_infection;
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

void GPU::ReporterUtils::output_genotype_frequency4(std::stringstream& ss, std::stringstream& ss2, const int& number_of_genotypes,
                                               GPU::PersonIndexByLocationStateAgeClass* pi) {
  auto sum1_all = 0.0;
  std::vector<double> result4_all(number_of_genotypes, 0.0);
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0; loc < number_of_locations; loc++) {
    std::vector<double> result4(number_of_genotypes, 0.0);
    std::vector<double> prmc4(number_of_genotypes, 0.0);
    std::vector<double> prmc4_all(number_of_genotypes, 0.0);
    auto sum1 = 0.0;

    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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
//            individual_genotype_cr[parasite_population->genotype()->genotype_id] = parasite_population->genotype()->daily_fitness_multiple_infection;
          }

          for (const auto genotype : individual_genotype_map) {
            result4[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
            result4_all[genotype.first] +=
                genotype.second / static_cast<double>(person->all_clonal_parasite_populations()->parasites()->size());
          }
        }
      }
    }
    // output per location
    // TODO: implement dynamic way to output for each location

//    for (auto& i : result3) {
//      i /= sum1;
//      ss << (sum1 == 0 ? 0 : i) << sep;
//    }

    std::map<int, int> prmc_genotype_map;
    auto tracking_day = Model::GPU_SCHEDULER->current_time() % Model::CONFIG->number_of_tracking_days();
    int sum_nullptr = 0;
    for(auto* prmc_genotypes : Model::GPU_MOSQUITO->genotypes_table[tracking_day][loc]) {
        if (prmc_genotypes == nullptr) {
            sum_nullptr++;
            continue;
        }
    }
    if (sum_nullptr == 0){
        for(auto* prmc_genotypes : Model::GPU_MOSQUITO->genotypes_table[tracking_day][loc]) {
          const auto g_id = prmc_genotypes->genotype_id;
          if (prmc_genotype_map.find(g_id) == prmc_genotype_map.end()) {
              prmc_genotype_map[prmc_genotypes->genotype_id] = 1;
          }
          else{
              prmc_genotype_map[prmc_genotypes->genotype_id] += 1;
          }
        }
        for (const auto genotype : prmc_genotype_map) {
            prmc4[genotype.first] += genotype.second / static_cast<double>(Model::GPU_MOSQUITO->genotypes_table[tracking_day][loc].size());
            prmc4_all[genotype.first] += genotype.second / static_cast<double>(Model::GPU_MOSQUITO->genotypes_table[tracking_day][loc].size());
        }

        for (auto& i : prmc4) {
            ss2 << i << sep;
        }
    }
  }
//  ss << group_sep;

  // this is for all locations
  for (auto& i : result4_all) {
    i /= sum1_all;
    ss << (sum1_all == 0 ? 0 : i) << sep;
  }
//
//  ss << group_sep;
//  ss << sum1_all << sep;
}

void GPU::ReporterUtils::output_3_genotype_frequency(std::stringstream& ss, const int& number_of_genotypes,
                                                GPU::PersonIndexByLocationStateAgeClass* pi) {
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

    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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

void GPU::ReporterUtils::output_moi(std::stringstream& ss, GPU::PersonIndexByLocationStateAgeClass* pi) {
  const auto number_of_locations = pi->vPerson().size();
  const auto number_of_age_classes = pi->vPerson()[0][0].size();

  for (auto loc = 0ul; loc < number_of_locations; loc++) {
    auto pop_sum_location = 0;
    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
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

void GPU::ReporterUtils::initialize_moi_file_logger() {
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
//    ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_treatment_by_location()[loc] << sep;
//  }
//
//  ss << group_sep;
//
//  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
//    ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] << sep;
//  }
//}