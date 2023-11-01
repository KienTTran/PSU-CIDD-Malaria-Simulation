/* 
 * File:   PkPdReporter.cpp
 * Author: Merlin
 * 
 * Created on October 29, 2014, 12:56 PM
 */

#include "AppInput.h"
#include "PkPdReporter.h"
#include "Model.h"
#include "Core/Random.h"
#include "MDC/ModelDataCollector.h"
#include "Core/Config/Config.h"
#include "Population/Population.h"
#include "Population/Properties/PersonIndexAll.h"
#include "Population/Person.h"
#include "Population/ClonalParasitePopulation.h"
#include "Population/SingleHostClonalParasitePopulations.h"
#include "Population/DrugsInBlood.h"
#include "Therapies/Drug.h"
#include "Population/ImmuneSystem.h"

PkPdReporter::PkPdReporter(AppInput* appInput) : appInput { appInput } {
  if (appInput && !appInput->output_file.empty()) {
    outputFStream.open(appInput->output_file);
  }
}

PkPdReporter::~PkPdReporter(){
  if (outputFStream.is_open()) {
    outputFStream.close();
  }
}

void PkPdReporter::initialize() {
}

void PkPdReporter::before_run() {
  //    std::cout << Model::RANDOM->seed() << std::endl;
}

void PkPdReporter::begin_time_step() {
  ss << Model::SCHEDULER->current_time();

  // parasiteamia profile

  for (int i = 0; i < Model::POPULATION->all_persons()->vPerson().size(); i++) {
    auto p_person = Model::POPULATION->all_persons()->vPerson()[i];
    if (p_person->all_clonal_parasite_populations()->size() > 0) {
      ss << sep << p_person->all_clonal_parasite_populations()->parasites()->at(0)->get_log10_infectious_density();
    } else {
      ss << sep << Model::CONFIG->parasite_density_level().log_parasite_density_cured;
    }
  }

//  auto start_drug_id = appInput->is_art ? 0 : 1;
//  // drug concentrate profile
//  for (int drug_id = start_drug_id; drug_id < appInput->number_of_drugs_in_combination + start_drug_id; drug_id++) {
//    for (int i = 0; i < Model::POPULATION->all_persons()->vPerson().size(); i++) {
//      auto p_person = Model::POPULATION->all_persons()->vPerson()[i];
//      bool found { false };
//      for (const auto& [did, active_drug] : *(p_person->drugs_in_blood()->drugs())) {
//        if (drug_id == did) {
//          found = true;
//          ss << sep << active_drug->last_update_value();
//          break;
//        }
//      }
//      if (!found) {
//        ss << sep << 0;
//      }
//    }
//  }
//
//  // drug killing profile
//  for (int drug_id = start_drug_id; drug_id < appInput->number_of_drugs_in_combination + start_drug_id; drug_id++) {
//    for (int i = 0; i < Model::POPULATION->all_persons()->vPerson().size(); i++) {
//      auto p_person = Model::POPULATION->all_persons()->vPerson()[i];
//
//      bool found { false };
//      for (const auto& [did, active_drug] : *(p_person->drugs_in_blood()->drugs())) {
//        if (drug_id == did) {
//          found = true;
//          ss << sep << active_drug->get_parasite_killing_rate(0);
//          break;
//        }
//      }
//      if (!found) {
//        ss << sep << 0;
//      }
//    }
//  }

  outputFStream << ss.str() << std::endl;
  ss.str("");
}

void PkPdReporter::after_time_step() {

  //    if (Model::SCHEDULER->current_time() % Model::CONFIG->report_frequency() == 0) {
  //
  //        std::cout << Model::SCHEDULER->current_time() << "\t";
  //        std::cout << Model::DATA_COLLECTOR->blood_slide_prevalence_by_location()[0]*100 << "\t";
  //        std::cout << Model::DATA_COLLECTOR->number_of_positive_by_location()[0]*100 / (double) Model::DATA_COLLECTOR->popsize_by_location()[0] << "\t";
  //
  //        accumulator_set<double, stats<tag::median, tag::mean, tag::variance > > acc;
  //        accumulator_set<double, stats<tag::tail_quantile<right> > > acc_q(tag::tail<right>::cache_size = 10000);
  //
  //        DoubleVector today_density;
  //        //        std::cout <<  Model::POPULATION->all_persons()->vPerson().size()<<std::endl;
  //        for (int i = 0; i < Model::POPULATION->all_persons()->vPerson().size(); i++) {
  //            Person* person = Model::POPULATION->all_persons()->vPerson()[i];
  //            if (person->all_clonal_parasite_populations()->parasites()->size() >= 1) {
  //                acc(person->all_clonal_parasite_populations()->parasites()->at(0)->last_update_log10_parasite_density());
  //                acc_q(person->all_clonal_parasite_populations()->parasites()->at(0)->last_update_log10_parasite_density());
  //                today_density.push_back(person->all_clonal_parasite_populations()->parasites()->at(0)->last_update_log10_parasite_density());
  //            } else {
  //                acc(Model::CONFIG->parasite_density_level().log_parasite_density_cured);
  //                acc_q(Model::CONFIG->parasite_density_level().log_parasite_density_cured);
  //                today_density.push_back(Model::CONFIG->parasite_density_level().log_parasite_density_cured);
  //            }
  //        }
  ////        std::cout << today_density.size() << std::endl;
  //        assert(today_density.size() == 10000);
  //
  //        accumulator_set<double, stats<tag::median, tag::mean, tag::variance > > acc_PRR;
  //        accumulator_set<double, stats<tag::tail_quantile<right> > > acc_q_PRR(tag::tail<right>::cache_size = 10000);
  //        DoubleVector PRR;
  //        for (int i = 0; i < today_density.size(); i++) {
  //            if (Model::SCHEDULER->current_time() == 0) {
  //                PRR.push_back(0);
  //                acc_PRR(0);
  //                acc_q_PRR(0);
  //            } else {
  //                double t_density = today_density[i];
  //                double y_density = yesterday_density_[i];
  //                if (today_density[i] > Model::CONFIG->parasite_density_level().log_parasite_density_cured) {
  //                    PRR.push_back(pow(10, y_density - t_density));
  //                    acc_PRR(pow(10, y_density - t_density));
  //                    acc_q_PRR(pow(10, y_density - t_density));
  //                }
  //            }
  //        }
  //
  //        yesterday_density_ = today_density;
  //
  //
  //        std::cout << mean(acc) << "\t";
  //        std::cout << sqrt(variance(acc)) << "\t";
  //
  //        std::cout << median(acc) << "\t";
  //        std::cout << quantile(acc_q, quantile_probability = 0.75) << "\t";
  //        std::cout << quantile(acc_q, quantile_probability = 0.25) << "\t";
  //
  //        std::cout << median(acc_PRR) << "\t";
  //        std::cout << quantile(acc_q_PRR, quantile_probability = 0.75) << "\t";
  //        std::cout << quantile(acc_q_PRR, quantile_probability = 0.25) << "\t";
  //
  //        std::cout << std::endl;
  //    }


}

void PkPdReporter::monthly_report() {}

void PkPdReporter::after_run() {
  Model::DATA_COLLECTOR->update_after_run();

  //    std::cout << 1-Model::DATA_COLLECTOR->current_TF_by_location()[0] << std::endl;
  // std::cout << Model::DATA_COLLECTOR->total_number_of_treatments_60_by_location()[0][1] << std::endl;
}
