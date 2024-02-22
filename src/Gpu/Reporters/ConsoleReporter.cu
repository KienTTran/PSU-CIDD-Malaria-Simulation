/* 
 * File:   ConsoleReporter.cu
 * Author: Merlin
 *
 * Created on August 1, 2013, 12:15 PM
 */

#include <iostream>
#include "fmt/printf.h"
#include "ConsoleReporter.cuh"
#include "Model.h"
#include "Core/Random.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Constants.h"
#include "ReporterUtils.cuh"

GPU::ConsoleReporter::ConsoleReporter() {
}

GPU::ConsoleReporter::~ConsoleReporter() {
}

void GPU::ConsoleReporter::initialize() {
}

void GPU::ConsoleReporter::before_run() {
  std::cout << "Seed:" << Model::RANDOM->seed() << std::endl;

}

void report_number_by_state(const int &location, GPU::PersonIndexByLocationStateAgeClass* pi) {
  //    std::cout << std::setw(10) << std::setprecision(3);
  for (int hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
    //        int sum = 0;
    //        for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
    //            sum += pi->vPerson()[location][hs][ac].size();
    //        }
    double v = Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[location][hs] * 100 /
               (double) Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];
    //        double v = sum;

    fmt::printf("%.3f\t", v);

  }

}

void GPU::ConsoleReporter::after_run() {
  std::cout << "==========================================================================" << std::endl;

  //total time
  double total_time_in_years = (Model::GPU_SCHEDULER->current_time() - Model::CONFIG->start_collect_data_day()) /
                               (double) Constants::DAYS_IN_YEAR();
  std::cout << "Total time (from equilibrium) : " << total_time_in_years << " years" << std::endl;

  //report EIR
  std::cout << "EIR by location:" << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    std::cout << Model::GPU_DATA_COLLECTOR->EIR_by_location()[location] << "\t";
  }
  std::cout << std::endl;

  //total number of bites
  std::cout << "Number of infectious bites:" << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    std::cout << Model::GPU_DATA_COLLECTOR->total_number_of_bites_by_location()[location] << "\t";
  }
  std::cout << std::endl;

  std::cout << "Number of clinical episodes:" << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    std::cout << Model::GPU_DATA_COLLECTOR->cumulative_clinical_episodes_by_location()[location] << "\t";
  }
  std::cout << std::endl;

  std::cout << "Percentage of bites on top 20% bitten" << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    std::cout << Model::GPU_DATA_COLLECTOR->percentage_bites_on_top_20_by_location()[location] * 100 << "%" << "\t";
  }
  std::cout << std::endl;

  std::cout << "NTF by location: " << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    double location_NTF = Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[location] * 100 /
                          (double) Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];
    location_NTF /= total_time_in_years;

    std::cout << location_NTF << "\t";
  }
  std::cout << std::endl;

  std::cout << "Number of mutations by location: " << std::endl;
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    std::cout << Model::GPU_DATA_COLLECTOR->cumulative_mutants_by_location()[location] << "\t";
  }
  std::cout << std::endl;

  for (int t_id = 0; t_id < Model::CONFIG->gpu_therapy_db().size(); t_id++) {

    int nTreaments = Model::GPU_DATA_COLLECTOR->number_of_treatments_with_therapy_ID()[t_id];
    int nSuccess = Model::GPU_DATA_COLLECTOR->number_of_treatments_success_with_therapy_ID()[t_id];
    int nFail = Model::GPU_DATA_COLLECTOR->number_of_treatments_fail_with_therapy_ID()[t_id];
    double pSuccess = (nTreaments == 0) ? 0 : nSuccess * 100.0 / nTreaments;

    std::cout << "Number of patients (with non-resistant parasite) treated with therapy " << t_id
              << " (% success) = "
              << nTreaments << " (" << pSuccess << "%)" << std::endl;
    std::cout << "Number Failed: " << nFail << "-" << nSuccess << "-" << nTreaments << std::endl;

  }

  std::cout << "Strategy UTL: " << Model::GPU_DATA_COLLECTOR->current_utl_duration() << std::endl;

  std::cout << "AMU per parasite population: " << Model::GPU_DATA_COLLECTOR->AMU_per_parasite_pop() << std::endl;
  std::cout << "AMU per per: " << Model::GPU_DATA_COLLECTOR->AMU_per_person() << std::endl;
  std::cout << "EAMU count only clinical caused parasite: "
            << Model::GPU_DATA_COLLECTOR->AMU_for_clinical_caused_parasite()
            << std::endl;

}

void GPU::ConsoleReporter::begin_time_step() {
}

void GPU::ConsoleReporter::monthly_report() {

  if (Model::GPU_SCHEDULER->current_time() % Model::CONFIG->report_frequency() == 0) {
//        Model::GPU_DATA_COLLECTOR->perform_population_statistic();

    std::cout << Model::GPU_SCHEDULER->current_time() << "\t";

    auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();

    for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
      std::cout << "||\t";
      report_number_by_state(location, pi);
      std::cout << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[location] * 100 << "\t";
      std::cout << Model::GPU_DATA_COLLECTOR->total_immune_by_location()[location] / Model::GPU_POPULATION->size(location)
                << "\t";
      std::cout << Model::GPU_DATA_COLLECTOR->current_RITF_by_location()[location] << "-"
                << Model::GPU_DATA_COLLECTOR->current_TF_by_location()[location] << "\t";
    }
    std::cout << std::endl;
  }
}
