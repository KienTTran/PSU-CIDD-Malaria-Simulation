//
// Created by Nguyen Tran on 3/5/2018.
//

#include "MMCReporter.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Core/Random.h"
#include "Constants.h"
#include "easylogging++.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Strategies/IStrategy.cuh"
#include "ReporterUtils.cuh"

GPU::MMCReporter::MMCReporter() = default;

void GPU::MMCReporter::initialize() {
  ReporterUtils::initialize_moi_file_logger();

}

void GPU::MMCReporter::before_run() {
  // // std::cout << "MMC Reporter" << std::endl;
  // for (auto genotype : (*Model::CONFIG->gpu_genotype_db)){
  //   std::cout << *genotype.second << std::endl;
  // }

}

void GPU::MMCReporter::begin_time_step() { }


void GPU::MMCReporter::print_treatment_failure_rate_by_therapy() {
  for (auto tf_by_therapy : Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()) {
    ss << tf_by_therapy << sep;
  }
}

void GPU::MMCReporter::print_ntf_by_location() {
  double sum_ntf = 0.0;
  ul pop_size = 0;
  for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    sum_ntf += Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[location];
    pop_size += Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];

  }

  ss << (sum_ntf * 100.0 / pop_size) << sep;
}

void GPU::MMCReporter::monthly_report() {
  ss << Model::GPU_SCHEDULER->current_time() << sep;
  ss << std::chrono::system_clock::to_time_t(Model::GPU_SCHEDULER->calendar_date) << sep;
  ss << date::format("%Y\t%m\t%d", Model::GPU_SCHEDULER->calendar_date) << sep;
  ss << Model::CONFIG->seasonal_info()->get_seasonal_factor(Model::GPU_SCHEDULER->calendar_date, 0) << sep;
  ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 1) << sep;
  ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 10) << sep;
  ss << Model::GPU_POPULATION->size() << sep;
  ss << group_sep;

  print_EIR_PfPR_by_location();
  ss << group_sep;
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_treatment_by_location()[loc] << sep;
  }
  ss << group_sep;
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] << sep;
  }
  ss << group_sep;

  ReporterUtils::output_genotype_frequency3(
      ss,
      Model::CONFIG->gpu_genotype_db.size(),
      Model::GPU_POPULATION->get_person_index<PersonIndexByLocationStateAgeClass>());

  ss << group_sep;
  print_ntf_by_location();
  ss << group_sep;
  print_treatment_failure_rate_by_therapy();
  ss << Model::GPU_DATA_COLLECTOR->current_TF_by_location()[0];
  CLOG(INFO, "monthly_reporter") << ss.str();
  ss.str("");
}

void GPU::MMCReporter::after_run() {
  ss.str("");
  ss << Model::RANDOM->seed() << sep << Model::CONFIG->number_of_locations() << sep;
  ss << Model::CONFIG->location_db()[0].beta << sep;
  ss << Model::CONFIG->location_db()[0].population_size << sep;
  print_EIR_PfPR_by_location();

  ss << group_sep;
  //output last strategy information
  ss << Model::GPU_TREATMENT_STRATEGY->id << sep;

  //output NTF
  const auto total_time_in_years = (Model::GPU_SCHEDULER->current_time() - Model::CONFIG->start_of_comparison_period()) /
                                   static_cast<double>(Constants::DAYS_IN_YEAR());

  auto sum_ntf = 0.0;
  ul pop_size = 0;
  for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    sum_ntf += Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[location];
    pop_size += Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];
  }

  ss << (sum_ntf * 100 / pop_size) / total_time_in_years << sep;
  ss << group_sep;

  ss << Model::GPU_DATA_COLLECTOR->cumulative_number_treatments_by_location()[0] << sep;
  ss << Model::GPU_DATA_COLLECTOR->cumulative_TF_by_location()[0] << sep;
  ss << Model::GPU_DATA_COLLECTOR->cumulative_clinical_episodes_by_location()[0] << sep;

  ss << group_sep;
  //print # mutation events of first 10 years
  int number_of_years = Model::GPU_DATA_COLLECTOR->number_of_mutation_events_by_year().size() >= 11 ? 11 :
                        Model::GPU_DATA_COLLECTOR->number_of_mutation_events_by_year().size();
  for (int i = 0; i < number_of_years; ++i) {
    ss << Model::GPU_DATA_COLLECTOR->number_of_mutation_events_by_year()[i] << sep;
  }

  CLOG(INFO, "summary_reporter") << ss.str();
  ss.str("");

  // Report MOI
  ReporterUtils::output_moi(ss, Model::GPU_POPULATION->get_person_index<PersonIndexByLocationStateAgeClass>());
}

void GPU::MMCReporter::print_EIR_PfPR_by_location() {
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
    //
    // EIR
    if (Model::GPU_DATA_COLLECTOR->EIR_by_location_year()[loc].empty()) {
      ss << 0 << sep;
    } else {
      ss << Model::GPU_DATA_COLLECTOR->EIR_by_location_year()[loc].back() << sep;
    }

    //pfpr <5 and all
    ss << Model::GPU_DATA_COLLECTOR->get_blood_slide_prevalence(loc, 0, 5) * 100 << sep;
    ss << Model::GPU_DATA_COLLECTOR->get_blood_slide_prevalence(loc, 2, 10) * 100 << sep;
    ss << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[loc] * 100 << sep;
//    std::cout << Model::GPU_POPULATION->size() << "\t"
//              << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[loc] * 100 << std::endl;
  }
}


