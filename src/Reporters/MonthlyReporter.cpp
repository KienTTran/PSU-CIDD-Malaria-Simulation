//
// Created by Nguyen Tran on 3/5/2018.
//

#include "MonthlyReporter.h"

#include <date/date.h>

#include "Constants.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Helpers/TimeHelpers.h"
#include "MDC/ModelDataCollector.h"
#include "Model.h"
#include "Population/Population.h"
#include "Population/Properties/PersonIndexByLocationStateAgeClass.h"
#include "ReporterUtils.h"
#include "Strategies/IStrategy.h"
#include "easylogging++.h"

MonthlyReporter::MonthlyReporter() = default;

MonthlyReporter::~MonthlyReporter() = default;

void MonthlyReporter::initialize() {
  gene_freq_file.open(fmt::format("{}/gene_freq_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
  monthly_data_file.open(fmt::format("{}/monthly_data_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
  summary_data_file.open(fmt::format("{}/summary_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
  gene_db_file.open(fmt::format("{}/gene_db_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
}

void MonthlyReporter::before_run() {}

void MonthlyReporter::begin_time_step() {}

void MonthlyReporter::monthly_report() {
    std::stringstream ss;

    ss << Model::SCHEDULER->current_time() << sep;
    ss << std::chrono::system_clock::to_time_t(Model::SCHEDULER->calendar_date) << sep;
    ss << date::format("%Y\t%m\t%d", Model::SCHEDULER->calendar_date) << sep;
    ss << Model::MODEL->get_seasonal_factor(Model::SCHEDULER->calendar_date, 0) << sep;
    ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 1) << sep;
    ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 10) << sep;
    ss << Model::POPULATION->size() << sep;
    ss << group_sep;

    print_EIR_PfPR_by_location(ss);
    ss << group_sep;
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    ss << Model::DATA_COLLECTOR->monthly_number_of_new_infections_by_location()[loc] << sep;
    }
    ss << group_sep;
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    ss << Model::DATA_COLLECTOR->monthly_number_of_treatment_by_location()[loc] << sep;
    }
    ss << group_sep;
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    ss << Model::DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] << sep;
    }
    ss << group_sep;
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::DATA_COLLECTOR->cumulative_NTF_by_location()[loc] << sep;
    }
    ss << group_sep;

  // including total number of positive individuals
  //  ReporterUtils::output_genotype_frequency3(ss, Model::CONFIG->genotype_db.size(),
  //                                            Model::POPULATION->get_person_index<PersonIndexByLocationStateAgeClass>());

  std::stringstream gene_freq_ss;
  ReporterUtils::output_genotype_frequency3(gene_freq_ss, Model::CONFIG->genotype_db.size(),
                                            Model::POPULATION->get_person_index<PersonIndexByLocationStateAgeClass>());

  gene_freq_file << gene_freq_ss.str() << std::endl;

  monthly_data_file << ss.str() << std::endl;
}

void MonthlyReporter::after_run() {
  std::stringstream ss;

  ss.str("");
  ss << Model::RANDOM->seed() << sep << Model::CONFIG->number_of_locations() << sep;
  ss << Model::CONFIG->location_db()[0].beta << sep;
  ss << Model::CONFIG->location_db()[0].population_size << sep;
  print_EIR_PfPR_by_location(ss);

  ss << group_sep;
  // output last strategy information
  ss << Model::TREATMENT_STRATEGY->id << sep;

  // output NTF
  const auto total_time_in_years = (Model::SCHEDULER->current_time() - Model::CONFIG->start_of_comparison_period())
                                   / static_cast<double>(Constants::DAYS_IN_YEAR());

  auto sum_ntf = 0.0;
  ul pop_size = 0;
  for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    sum_ntf += Model::DATA_COLLECTOR->cumulative_NTF_by_location()[location];
    pop_size += Model::DATA_COLLECTOR->popsize_by_location()[location];
  }
  ss << (sum_ntf * 100 / pop_size) / total_time_in_years << sep;

  ss << group_sep;
  for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++)
  {
    for (auto age = 0; age < 60; age++){
      ss << Model::DATA_COLLECTOR->cumulative_clinical_episodes_by_location_age()[location][age]/total_time_in_years/Model::DATA_COLLECTOR->popsize_by_location_age()[location][age] << sep;
    }
  }
  summary_data_file << ss.str() << std::endl;

  for (auto [g_id, genotype] : Model::CONFIG->genotype_db) {
    gene_db_file << g_id << sep << genotype->aa_sequence << std::endl;
  }

  gene_freq_file.close();
  monthly_data_file.close();
  summary_data_file.close();
  gene_db_file.close();
}

void MonthlyReporter::print_EIR_PfPR_by_location(std::stringstream& ss) {
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
    //
    // EIR
    if (Model::DATA_COLLECTOR->EIR_by_location_year()[loc].empty()) {
      ss << 0 << sep;
    } else {
      ss << Model::DATA_COLLECTOR->EIR_by_location_year()[loc].back() << sep;
    }
    ss << group_sep;
    // pfpr <5 , 2-10 and all
    ss << Model::DATA_COLLECTOR->get_blood_slide_prevalence(loc, 2, 10) * 100 << sep;
    ss << Model::DATA_COLLECTOR->get_blood_slide_prevalence(loc, 0, 5) * 100 << sep;
    ss << Model::DATA_COLLECTOR->blood_slide_prevalence_by_location()[loc] * 100 << sep;
  }
}
