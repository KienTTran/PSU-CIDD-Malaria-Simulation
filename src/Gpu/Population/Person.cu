/*
 * File:   Person.cpp
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:25 PM
 */

#include "Person.cuh"

#include <algorithm>
#include <cmath>
#include <glm/ext/matrix_transform.hpp>

#include "ClonalParasitePopulation.cuh"
#include "Constants.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Events/BirthdayEvent.cuh"
#include "Gpu/Events/CirculateToTargetLocationNextDayEvent.cuh"
#include "Gpu/Events/EndClinicalByNoTreatmentEvent.cuh"
#include "Gpu/Events/EndClinicalDueToDrugResistanceEvent.cuh"
#include "Gpu/Events/EndClinicalEvent.cuh"
#include "Gpu/Events/MatureGametocyteEvent.cuh"
#include "Gpu/Events/MoveParasiteToBloodEvent.cuh"
#include "Gpu/Events/ProgressToClinicalEvent.cuh"
#include "Gpu/Events/ReceiveTherapyEvent.cuh"
#include "Gpu/Events/ReturnToResidenceEvent.cuh"
#include "Gpu/Events/TestTreatmentFailureEvent.cuh"
#include "Gpu/Events/UpdateWhenDrugIsPresentEvent.cuh"
#include "Helpers/ObjectHelpers.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Model.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Therapies/DrugDatabase.cuh"
#include "Gpu/Therapies/Drug.cuh"
#include "Gpu/Therapies/DrugType.cuh"
#include "Gpu/Therapies/MACTherapy.cuh"
#include "Gpu/Therapies/SCTherapy.cuh"
#include "Helpers/UniqueId.hxx"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/DrugsInBlood.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"
#include "Gpu/Population/Population.cuh"
#include "ImmuneSystem.cuh"

GPU::Person::Person()
    : location_(-1),
      residence_location_(-1),
      host_state_(SUSCEPTIBLE),
      age_(-1),
      age_class_(-1),
      birthday_(-1),
      latest_update_time_(-1),
      innate_relative_biting_rate { 0 },
      moving_level_(-1),
      liver_parasite_type_(nullptr),
      number_of_times_bitten_(0),
      number_of_trips_taken_(0),
      last_therapy_id_(0) {
  population_ = nullptr;
  immune_system_ = nullptr;
  all_clonal_parasite_populations_ = nullptr;
  drugs_in_blood_ = nullptr;

  today_infections_ = nullptr;
  today_target_locations_ = nullptr;
  latest_update_time_ = -1;
  location_col_ = -1;
  location_row_ = -1;
  model_ = glm::mat4(1.0f);
  color_ = glm::vec4(0.0f);
  id_ = -1;
  index_ = -1;
}

void GPU::Person::init() {
  Dispatcher::init();

  immune_system_ = new GPU::ImmuneSystem(this);

  all_clonal_parasite_populations_ = new GPU::SingleHostClonalParasitePopulations(this);
  all_clonal_parasite_populations_->init();

  drugs_in_blood_ = new GPU::DrugsInBlood(this);
  drugs_in_blood_->init();

  today_infections_ = new IntVector();
  today_target_locations_ = new IntVector();

  person_index_gpu = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
}

GPU::Person::~Person() {
  Dispatcher::clear_events();
  ObjectHelpers::delete_pointer<GPU::ImmuneSystem>(immune_system_);
  ObjectHelpers::delete_pointer<GPU::SingleHostClonalParasitePopulations>(all_clonal_parasite_populations_);
  ObjectHelpers::delete_pointer<GPU::DrugsInBlood>(drugs_in_blood_);
  ObjectHelpers::delete_pointer<IntVector>(today_infections_);
  ObjectHelpers::delete_pointer<IntVector>(today_target_locations_);
}

void GPU::Person::NotifyChange(const GPU::Person::Property& property, const void* oldValue, const void* newValue) {
  if (population_ != nullptr) {
    population_->notify_change(this, property, oldValue, newValue);
  }
}

int GPU::Person::location() const {
  return location_;
}

void GPU::Person::set_location(const int& value) {
  if (location_ != value) {
    if (Model::GPU_DATA_COLLECTOR != nullptr) {
      const auto day_diff = (Constants::DAYS_IN_YEAR() - Model::GPU_SCHEDULER->current_day_in_year());
      if (location_ != -1) {
        Model::GPU_DATA_COLLECTOR->update_person_days_by_years(location_, -day_diff);
      }
      Model::GPU_DATA_COLLECTOR->update_person_days_by_years(value, day_diff);
    }

    Model::GPU_DATA_COLLECTOR->record_1_migration(this, location_, value);

    NotifyChange(LOCATION, &location_, &value);

    location_ = value;
  }
}


int GPU::Person::index() const {
  return index_;
}

void GPU::Person::set_index(const int& value) {
  index_ = value;
}

long GPU::Person::id() const {
    return id_;
}

void GPU::Person::set_id(const long& value) {
  id_ = value;
}

GPU::Person::HostStates GPU::Person::host_state() const {
  return host_state_;
}

void GPU::Person::set_host_state(const GPU::Person::HostStates& value) {
  if (host_state_ != value) {
    NotifyChange(HOST_STATE, &host_state_, &value);
    if (value == DEAD) {
      // clear also remove all infection forces
      all_clonal_parasite_populations_->clear();
      clear_events();

      //
      //            Model::STATISTIC->update_person_days_by_years(location_, -(Constants::DAYS_IN_YEAR() -
      //            Model::GPU_SCHEDULER->current_day_in_year()));
      Model::GPU_DATA_COLLECTOR->record_1_death(location_, birthday_, number_of_times_bitten_, age_class_, age_);
    }

    host_state_ = value;
  }
}

int GPU::Person::age() const {
  return age_;
}

void GPU::Person::set_age(const int& value) {
  if (age_ != value) {
    // TODO::if age access the limit of age structure i.e. 100, remove person???

    NotifyChange(AGE, &age_, &value);

    // update bitting rate level
    age_ = value;

    // update age class
    if (Model::MODEL != nullptr) {
      auto ac = age_class_ == -1 ? 0 : age_class_;

      while (ac < (Model::CONFIG->number_of_age_classes() - 1) && age_ >= Model::CONFIG->age_structure()[ac]) {
        ac++;
      }

      set_age_class(ac);
    }
  }
}

int GPU::Person::age_class() const {
  return age_class_;
}

void GPU::Person::set_age_class(const int& value) {
  if (age_class_ != value) {
    NotifyChange(AGE_CLASS, &age_class_, &value);
    age_class_ = value;
  }
}

int GPU::Person::moving_level() const {
  return moving_level_;
}

void GPU::Person::set_moving_level(const int& value) {
  if (moving_level_ != value) {
    NotifyChange(MOVING_LEVEL, &moving_level_, &value);
    moving_level_ = value;
  }
}

void GPU::Person::increase_age_by_1_year() {
  set_age(age_ + 1);
}

GPU::ImmuneSystem* GPU::Person::immune_system() const {
  return immune_system_;
}

void GPU::Person::set_immune_system(GPU::ImmuneSystem* value) {
  if (immune_system_ != value) {
    delete immune_system_;
    immune_system_ = value;
  }
}

void GPU::Person::set_birthday(const int &value) {
    birthday_ = value;
}

int GPU::Person::birthday() const {
    return birthday_;
}

GPU::ClonalParasitePopulation* GPU::Person::add_new_parasite_to_blood(GPU::Genotype* parasite_type) {
  auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
  auto* blood_parasite = new GPU::ClonalParasitePopulation(parasite_type);

  all_clonal_parasite_populations_->add(blood_parasite);
  /*
   * Need to set person and parasite index after adding parasite
   * Also the genotype is set here for the first time
   * */
  blood_parasite->set_person(this);
  blood_parasite->set_index(all_clonal_parasite_populations_->size() - 1);
  std::copy( parasite_type->aa_sequence.begin(),
             parasite_type->aa_sequence.end(),
             person_index_gpu->h_person_update_info()[index_].parasite_genotype[blood_parasite->index()]);
  person_index_gpu->h_person_update_info()[index_].parasite_genotype[blood_parasite->index()][MAX_GENOTYPE_LOCUS] = '\0';
  blood_parasite->set_index(all_clonal_parasite_populations_->size() - 1);
  person_index_gpu->h_person_update_info()[index_].parasite_genotype_fitness_multiple_infection[blood_parasite->index()]
  = parasite_type->daily_fitness_multiple_infection;

  blood_parasite->set_last_update_log10_parasite_density(Model::CONFIG->parasite_density_level().log_parasite_density_from_liver);

  LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
    << fmt::format("{} GPU::Person::add_new_parasite_to_blood before {} {} {}",
               index_,all_clonal_parasite_populations_->size()-1,
               parasite_type->aa_sequence.c_str(),
               blood_parasite->last_update_log10_parasite_density());

//  printf("add_new_parasite_to_blood: id %d <--> id %d p_index %d <--> c_index %d\n", id_,blood_parasite->id(),index_,blood_parasite->index());
  person_index_gpu->h_person_update_info()[index_].person_id = id_;
  person_index_gpu->h_person_update_info()[index_].person_index = index_;
  person_index_gpu->h_person_update_info()[index_].person_latest_update_time = latest_update_time_;
  person_index_gpu->h_person_update_info()[index_].parasites_size = all_clonal_parasite_populations_->size();
  person_index_gpu->h_person_update_info()[index_].parasites_current_index = blood_parasite->index();
  person_index_gpu->h_person_update_info()[index_].parasite_id[blood_parasite->index()] = blood_parasite->id();
  LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
    << fmt::format("{} GPU::Person::add_new_parasite_to_blood after {} {} {}",
               index_,all_clonal_parasite_populations_->size(),
               person_index_gpu->h_person_update_info()[index_].parasite_genotype[blood_parasite->index()],
               person_index_gpu->h_person_update_info()[index_].parasite_last_update_log10_parasite_density[blood_parasite->index()]);

  return blood_parasite;
}

double GPU::Person::relative_infectivity(const double& log10_parasite_density) {
  if (log10_parasite_density == ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY) return 0.0;

  // this sigma has already taken 'ln' and 'log10' into account
  const auto d_n = log10_parasite_density * Model::CONFIG->relative_infectivity().sigma
                   + Model::CONFIG->relative_infectivity().ro_star;
  const auto p = 0.001; //Model::RANDOM->cdf_standard_normal_distribution(d_n);

  return p * p + 0.01;
}

double GPU::Person::get_probability_progress_to_clinical() {
  return immune_system_->get_clinical_progression_probability(Model::CONFIG->immune_system_information(),latest_update_time_,Model::GPU_SCHEDULER->current_time());
}

void GPU::Person::cancel_all_other_progress_to_clinical_events_except(GPU::Event* event) const {
  for (auto* e : *events()) {
    if (e != event && dynamic_cast<GPU::ProgressToClinicalEvent*>(e) != nullptr) {
      //            std::cout << "Hello"<< std::endl;
      e->executable = false;
    }
  }
}

void GPU::Person::cancel_all_events_except(GPU::Event* event) const {
  for (auto* e : *events()) {
    if (e != event) {
      //            e->set_dispatcher(nullptr);
      e->executable = false;
    }
  }
}

// void GPU::Person::record_treatment_failure_for_test_treatment_failure_events() {
//
//     for(GPU::Event* e :  *events()) {
//         if (dynamic_cast<GPU::TestTreatmentFailureEvent*> (e) != nullptr && e->executable()) {
//             //            e->set_dispatcher(nullptr);
//             //record treatment failure
//             Model::GPU_DATA_COLLECTOR->record_1_treatment_failure_by_therapy(location_, age_,
//             ((GPU::TestTreatmentFailureEvent*) e)->therapyId());
//
//         }
//     }
// }

void GPU::Person::change_all_parasite_update_function(GPU::ParasiteDensityUpdateFunction* from,
                                                 GPU::ParasiteDensityUpdateFunction* to) {
  all_clonal_parasite_populations_->change_all_parasite_update_function(from, to);
}

bool GPU::Person::will_progress_to_death_when_receive_no_treatment() {
  // yes == death
  const auto p = Model::RANDOM->random_flat(0.0, 1.0);
  return p <= Model::CONFIG->mortality_when_treatment_fail_by_age_class()[age_class_];
}

bool GPU::Person::will_progress_to_death_when_recieve_treatment() {
  LOG_IF(age_class_ > 100,INFO) << id_ << " " << index_ << "Age is greater than 100";
  // yes == death
  double P = Model::RANDOM->random_flat(0.0, 1.0);
  // 90% lower than no treatment
  return P <= Model::CONFIG->mortality_when_treatment_fail_by_age_class()[age_class_] * (1 - 0.9);
}

void GPU::Person::schedule_progress_to_clinical_event_by(GPU::ClonalParasitePopulation* blood_parasite) {
  const auto time =
      (age_ <= 5) ? Model::CONFIG->days_to_clinical_under_five() : Model::CONFIG->days_to_clinical_over_five();

  GPU::ProgressToClinicalEvent::schedule_event(Model::GPU_SCHEDULER, this, blood_parasite,
                                          Model::GPU_SCHEDULER->current_time() + time);
}

void GPU::Person::schedule_test_treatment_failure_event(GPU::ClonalParasitePopulation* blood_parasite, const int& testing_day,
                                                   const int& t_id) {
  GPU::TestTreatmentFailureEvent::schedule_event(Model::GPU_SCHEDULER, this, blood_parasite,
                                            Model::GPU_SCHEDULER->current_time() + testing_day, t_id);
}

int GPU::Person::complied_dosing_days(const int& dosing_day) const {
  if (Model::CONFIG->p_compliance() < 1) {
    const auto p = Model::RANDOM->random_flat(0.0, 1.0);
    if (p > Model::CONFIG->p_compliance()) {
      // do not comply
      const auto a = (Model::CONFIG->min_dosing_days() - dosing_day) / (1 - Model::CONFIG->p_compliance());
      return static_cast<int>(std::ceil(a * p + Model::CONFIG->min_dosing_days() - a));
    }
  }
  return dosing_day;
}

void GPU::Person::receive_therapy(GPU::Therapy* therapy, GPU::ClonalParasitePopulation* clinical_caused_parasite,
                             bool is_part_of_MAC_therapy) {
  // if therapy is SCTherapy
  auto* sc_therapy = dynamic_cast<GPU::SCTherapy*>(therapy);
  if (sc_therapy != nullptr) {
    for (int j = 0; j < sc_therapy->drug_ids.size(); ++j) {
      auto dosing_days = sc_therapy->drug_ids.size() == sc_therapy->dosing_day.size() ? sc_therapy->dosing_day[j]
                                                                                      : sc_therapy->dosing_day[0];

      dosing_days = complied_dosing_days(dosing_days);
      int drug_id = sc_therapy->drug_ids[j];
      LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
        << fmt::format("{} receive_therapy: SCTherapy\n",index_);
      add_drug_to_blood(Model::CONFIG->gpu_drug_db()->at(drug_id), dosing_days, is_part_of_MAC_therapy);
    }
  } else {
    // else if therapy is MACTherapy
    auto* mac_therapy = dynamic_cast<GPU::MACTherapy*>(therapy);
    starting_drug_values_for_MAC.clear();
    assert(mac_therapy != nullptr);
    for (auto i = 0; i < mac_therapy->therapy_ids().size(); i++) {
      const auto therapy_id = mac_therapy->therapy_ids()[i];
      const auto start_day = mac_therapy->start_at_days()[i];

      if (start_day == 1) {
        LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
          << fmt::format("{} receive_therapy: MACTherapy\n",index_);
        receive_therapy(Model::CONFIG->gpu_therapy_db()[therapy_id], clinical_caused_parasite, true);
      } else {
        assert(start_day > 1);
        GPU::ReceiveTherapyEvent::schedule_event(Model::GPU_SCHEDULER, this, Model::CONFIG->gpu_therapy_db()[therapy_id],
                                            Model::GPU_SCHEDULER->current_time() + start_day - 1, clinical_caused_parasite,
                                            true);
      }
    }
  }
  last_therapy_id_ = therapy->id();
}

void GPU::Person::add_drug_to_blood(GPU::DrugType* dt, const int& dosing_days, bool is_part_of_MAC_therapy) {
  auto* drug = new GPU::Drug(dt);
  drug->set_dosing_days(dosing_days);
  drug->set_last_update_time(Model::GPU_SCHEDULER->current_time());

  const auto sd = dt->age_group_specific_drug_concentration_sd()[age_class_];
  const auto mean_drug_absorption = dt->age_specific_drug_absorption()[age_class_];
  double drug_level = Model::RANDOM->random_normal_truncated(mean_drug_absorption, sd);

  if (is_part_of_MAC_therapy) {
    if (drugs_in_blood()->drugs()->find(dt->id()) != drugs_in_blood()->drugs()->end()) {
      // long-half life drugs
      drug_level = drugs_in_blood()->get_drug(dt->id())->starting_value();
    } else if (starting_drug_values_for_MAC.find(dt->id()) != starting_drug_values_for_MAC.end()) {
      // short half-life drugs
      drug_level = starting_drug_values_for_MAC[dt->id()];
    }
    // store the override value or the default one
    starting_drug_values_for_MAC[dt->id()] = drug_level;
  }
  drug->set_starting_value(drug_level);

  if (drugs_in_blood_->is_drug_in_blood(dt)) {
    drug->set_last_update_value(drugs_in_blood_->get_drug(dt->id())->last_update_value());
  } else {
    drug->set_last_update_value(0.0);
  }

  drug->set_start_time(Model::GPU_SCHEDULER->current_time());
  drug->set_end_time(Model::GPU_SCHEDULER->current_time() + dt->get_total_duration_of_drug_activity(dosing_days));
  drugs_in_blood_->add_drug(drug);
}

void GPU::Person::schedule_update_by_drug_event(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  GPU::UpdateWhenDrugIsPresentEvent::schedule_event(Model::GPU_SCHEDULER, this, clinical_caused_parasite,
                                               Model::GPU_SCHEDULER->current_time() + 1);
}

void GPU::Person::schedule_end_clinical_event(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  int dClinical = Model::RANDOM->random_normal(7, 2);
  dClinical = std::min<int>(std::max<int>(dClinical, 5), 14);

  GPU::EndClinicalEvent::schedule_event(Model::GPU_SCHEDULER, this, clinical_caused_parasite,
                                   Model::GPU_SCHEDULER->current_time() + dClinical);
}

void GPU::Person::schedule_end_clinical_by_no_treatment_event(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  auto d_clinical = Model::RANDOM->random_normal(7, 2);
  d_clinical = std::min<int>(std::max<int>(d_clinical, 5), 14);

  GPU::EndClinicalByNoTreatmentEvent::schedule_event(Model::GPU_SCHEDULER, this, clinical_caused_parasite,
                                                Model::GPU_SCHEDULER->current_time() + d_clinical);
}

void GPU::Person::change_state_when_no_parasite_in_blood() {
  if (all_clonal_parasite_populations_->size() == 0) {
    if (liver_parasite_type_ == nullptr) {
      //        std::cout << "S" << std::endl;
      //        std::cout << host_state_<< std::endl;
      set_host_state(SUSCEPTIBLE);
      //        std::cout << "ES" << std::endl;

    } else {
      set_host_state(EXPOSED);
    }
    immune_system_->set_increase(false);
  }
}

void GPU::Person::determine_relapse_or_not(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  if (all_clonal_parasite_populations_->contain(clinical_caused_parasite)) {
    const auto p = Model::RANDOM->random_flat(0.0, 1.0);

    if (p <= Model::CONFIG->p_relapse()) {
      //        if (P <= get_probability_progress_to_clinical()) {
      // progress to clinical after several days
      clinical_caused_parasite->set_update_function(Model::MODEL->gpu_progress_to_clinical_update_function());
      clinical_caused_parasite->set_last_update_log10_parasite_density(
          Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic);
      schedule_relapse_event(clinical_caused_parasite, Model::CONFIG->relapse_duration());

    } else {
      // progress to clearance
      if (clinical_caused_parasite->last_update_log10_parasite_density()
          > Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic) {
        clinical_caused_parasite->set_last_update_log10_parasite_density(
            Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic);
      }
      clinical_caused_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
    }
  }
}

void GPU::Person::determine_clinical_or_not(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  if (all_clonal_parasite_populations_->contain(clinical_caused_parasite)) {
    const auto p = Model::RANDOM->random_flat(0.0, 1.0);

    //        if (P <= Model::CONFIG->p_relapse()) {
    //    if (Model::GPU_SCHEDULER->current_time() >= 2000 && Model::GPU_SCHEDULER->current_time() <= 2010)
    //      std::cout << this->age() << "\t" << this->immune_system()->get_current_value() << "\t"
    //                << get_probability_progress_to_clinical()
    //                << std::endl;
    if (p <= get_probability_progress_to_clinical()) {
      // progress to clinical after several days
      clinical_caused_parasite->set_update_function(Model::MODEL->gpu_progress_to_clinical_update_function());
      clinical_caused_parasite->set_last_update_log10_parasite_density(
          Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic);
      schedule_relapse_event(clinical_caused_parasite, Model::CONFIG->relapse_duration());

    } else {
      // progress to clearance

      clinical_caused_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
    }
  }
}

void GPU::Person::schedule_relapse_event(GPU::ClonalParasitePopulation* clinical_caused_parasite, const int& time_until_relapse) {
  int duration = Model::RANDOM->random_normal(time_until_relapse, 15);
  duration = std::min<int>(std::max<int>(duration, time_until_relapse - 15), time_until_relapse + 15);
  GPU::ProgressToClinicalEvent::schedule_event(Model::GPU_SCHEDULER, this, clinical_caused_parasite,
                                          Model::GPU_SCHEDULER->current_time() + duration);
}

void GPU::Person::update() {
  //    std::cout << "Person Update"<< std::endl;
  // already update
  assert(host_state_ != DEAD);

  if (latest_update_time_ == Model::GPU_SCHEDULER->current_time()) return;

//  for(auto *parasite: *all_clonal_parasite_populations_->parasites()){
//    LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//      << fmt::format("{} CPU update_all_individuals before update parasite {} {} {} {} {} {}",
//             index_,
//             parasite->index(),
//             all_clonal_parasite_populations_->latest_update_time(),
//             Model::GPU_SCHEDULER->current_time(),
//             parasite->update_function()->type(),
//             parasite->genotype()->aa_sequence.c_str(),
//             parasite->last_update_log10_parasite_density());
//  }
  //    std::cout << "ppu"<< std::endl;
  // update the density of each blood parasite in parasite population
  // parasite will be killed by immune system

  all_clonal_parasite_populations_->update();

//  for (auto &drug : *drugs_in_blood_->drugs()) {
//    LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//      << fmt::format("{} CPU update_all_individuals before update drug {} {} {} {} {}",
//             index_,drug.second->start_time(),drug.second->last_update_time(),drug.first,
//             drug.second->starting_value(),drug.second->last_update_value());
//  }

//  printf("DrugsInBlood::update drugs_in_blood_ size %d\n", drugs_in_blood_->size());
  // update all drugs concentration
  drugs_in_blood_->update();

//  for (auto &drug : *drugs_in_blood_->drugs()) {
//    LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//      << fmt::format("{} CPU update_all_individuals after update drug {} {} {} {} {}",
//             index_,drug.second->start_time(),drug.second->last_update_time(),drug.first,
//             drug.second->starting_value(),drug.second->last_update_value());
//  }

  // update drug activity on parasite
  all_clonal_parasite_populations_->update_by_drugs(drugs_in_blood_);

//  LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//    << fmt::format("{} CPU update_all_individuals before update immune {} {} {}",
//                   index_,
//                   latest_update_time(),
//                   Model::GPU_SCHEDULER->current_time(),
//                   immune_system()->get_lastest_immune_value());

  immune_system_->update(Model::CONFIG->immune_system_information(),latest_update_time_,Model::GPU_SCHEDULER->current_time());

  update_current_state();

  // update bitting level only less than 1 to save performance
  //  the other will be update in birthday event
  update_relative_bitting_rate();

  latest_update_time_ = Model::GPU_SCHEDULER->current_time();
  //    std::cout << "End Person Update"<< std::endl;

//  for(auto *parasite: *all_clonal_parasite_populations_->parasites()){
//      LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//        << fmt::format("{} CPU update_all_individuals after update parasite {} {} {} {} {} {}",
//             index_,
//             parasite->index(),
//             all_clonal_parasite_populations_->latest_update_time(),
//             Model::GPU_SCHEDULER->current_time(),
//             parasite->update_function()->type(),
//             parasite->genotype()->aa_sequence.c_str(),
//             parasite->last_update_log10_parasite_density());
//  }

//  for (auto &drug : *drugs_in_blood_->drugs()) {
//    LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//      << fmt::format("{} CPU update_all_individuals after update drug clear {} {} {} {} {}",
//             index_,drug.second->start_time(),drug.second->last_update_time(),drug.first,
//             drug.second->starting_value(),drug.second->last_update_value());
//  }

//  LOG_IF(index_ >= 1040 && index_ <= 1045,INFO)
//    << fmt::format("{} CPU update_all_individuals after update immune {} {} {}",
//                   index_,
//                   latest_update_time(),
//                   Model::GPU_SCHEDULER->current_time(),
//                   immune_system()->get_lastest_immune_value());
}

void GPU::Person::update_relative_bitting_rate() {
  if (Model::CONFIG->using_age_dependent_bitting_level()) {
    current_relative_biting_rate = innate_relative_biting_rate * get_age_dependent_biting_factor();
  } else {
    current_relative_biting_rate = innate_relative_biting_rate;
  }
}

void GPU::Person::update_current_state() {
  //    std::cout << "ccod" << std::endl;
  // clear drugs <=0.1
  drugs_in_blood_->clear_cut_off_drugs_by_event(nullptr);
  //    std::cout << "ccp" << std::endl;
  // clear cured parasite
  all_clonal_parasite_populations_->clear_cured_parasites();

  //    std::cout << "change state" << std::endl;
  if (all_clonal_parasite_populations_->size() == 0) {
    change_state_when_no_parasite_in_blood();
  } else {
    immune_system_->set_increase(true);
  }
}

void GPU::Person::randomly_choose_parasite() {
  if (today_infections_->empty()) {
    // already chose
    return;
  }
  if (today_infections_->size() == 1) {
    infected_by(today_infections_->at(0));
  } else {
    const std::size_t index_random_parasite = Model::RANDOM->random_uniform(today_infections_->size());
    infected_by(today_infections_->at(index_random_parasite));
  }

  today_infections_->clear();
}

void GPU::Person::infected_by(const int& parasite_type_id) {
  // only infect if liver is available :D
  if (liver_parasite_type_ == nullptr) {
    if (host_state_ == SUSCEPTIBLE) {
      set_host_state(EXPOSED);
    }

    GPU::Genotype* genotype = Model::CONFIG->gpu_genotype_db.at(parasite_type_id);
    set_liver_parasite_type(genotype);
    std::copy( genotype->aa_sequence.begin(),
               genotype->aa_sequence.end(),
               person_index_gpu->h_person_update_info()[index_].person_liver_parasite_genotype);
    person_index_gpu->h_person_update_info()[index_].person_liver_parasite_genotype[MAX_GENOTYPE_LOCUS] = '\0';

    // move parasite to blood in next 7 days
    schedule_move_parasite_to_blood(genotype, 7);
  }
}

void GPU::Person::schedule_move_parasite_to_blood(GPU::Genotype* genotype, const int& time) {
  GPU::MoveParasiteToBloodEvent::schedule_event(Model::GPU_SCHEDULER, this, genotype, Model::GPU_SCHEDULER->current_time() + time);
}

void GPU::Person::schedule_mature_gametocyte_event(GPU::ClonalParasitePopulation* clinical_caused_parasite) {
  const auto day_mature_gametocyte = (age_ <= 5) ? Model::CONFIG->days_mature_gametocyte_under_five()
                                                 : Model::CONFIG->days_mature_gametocyte_over_five();
  GPU::MatureGametocyteEvent::schedule_event(Model::GPU_SCHEDULER, this, clinical_caused_parasite,
                                        Model::GPU_SCHEDULER->current_time() + day_mature_gametocyte);
}

void GPU::Person::randomly_choose_target_location() {
  if (today_target_locations_->empty()) {
    // already chose
    return;
  }

  auto target_location { -1 };
  if (today_target_locations_->size() == 1) {
    target_location = today_target_locations_->at(0);
  } else {
    const int index_random_location = Model::RANDOM->random_uniform(today_target_locations_->size());
    target_location = today_target_locations_->at(index_random_location);
  }

  schedule_move_to_target_location_next_day_event(target_location);

  today_target_locations_->clear();
}

void GPU::Person::schedule_move_to_target_location_next_day_event(const int& location) {
  this->number_of_trips_taken_ += 1;
  GPU::CirculateToTargetLocationNextDayEvent::schedule_event(Model::GPU_SCHEDULER, this, location,
                                                        Model::GPU_SCHEDULER->current_time() + 1);
}

bool GPU::Person::has_return_to_residence_event() const {
  for (GPU::Event* e : *events()) {
    if (dynamic_cast<GPU::ReturnToResidenceEvent*>(e) != nullptr) {
      return true;
    }
  }
  return false;
}

void GPU::Person::cancel_all_return_to_residence_events() const {
  for (GPU::Event* e : *events()) {
    if (dynamic_cast<GPU::ReturnToResidenceEvent*>(e) != nullptr) {
      e->executable = false;
    }
  }
}

bool GPU::Person::has_detectable_parasite() const {
  return all_clonal_parasite_populations_->has_detectable_parasite();
}

void GPU::Person::increase_number_of_times_bitten() {
  if (Model::GPU_SCHEDULER->current_time() >= Model::CONFIG->start_collect_data_day()) {
    number_of_times_bitten_++;
  }
}

void GPU::Person::move_to_population(GPU::Population* target_population) {
  assert(population_ != target_population);

  population_->remove_person(this);
  target_population->add_person(this);
}

bool GPU::Person::has_birthday_event() const {
  for (GPU::Event* e : *events()) {
    if (dynamic_cast<GPU::BirthdayEvent*>(e) != nullptr) {
      return true;
    }
  }
  return false;
}

bool GPU::Person::has_update_by_having_drug_event() const {
  for (GPU::Event* e : *events()) {
    if (dynamic_cast<GPU::UpdateWhenDrugIsPresentEvent*>(e) != nullptr) {
      return true;
    }
  }
  return false;
}

double GPU::Person::get_age_dependent_biting_factor() const {
  //
  // 0.00 - 0.25  -  6.5
  // 0.25 - 0.50  -  8.0
  // 0.50 - 0.75  -  9.0
  // 0.75 - 1.00  -  9.5
  // 1.00 - 2.00  -  11.0
  // 2.00 - 3.00  -  13.5
  // 3.00 - 4.00  -  15.5
  // 4.00 - 5.00  -  17.5
  // + 2.75kg until 20
  // then divide by 61.5

  if (age_ < 1) {
    const auto age = ((Model::GPU_SCHEDULER->current_time() - birthday_) % Constants::DAYS_IN_YEAR())
                     / static_cast<double>(Constants::DAYS_IN_YEAR());
    if (age < 0.25) return 0.106;
    if (age < 0.5) return 0.13;
    if (age < 0.75) return 0.1463;
    return 0.1545;
  }
  if (age_ < 2) return 0.1789;
  if (age_ < 3) return 0.2195;
  if (age_ < 4) return 0.2520;
  if (age_ < 20) return (17.5 + (age_ - 4) * 2.75) / 61.5;
  return 1.0;
}

double GPU::Person::p_infection_from_an_infectious_bite() const {
  return (1 - immune_system_->get_current_value(Model::CONFIG->immune_system_information(),latest_update_time_,Model::GPU_SCHEDULER->current_time())) / 8.333 + 0.04;
}

bool GPU::Person::isGametocytaemic() const {
  return all_clonal_parasite_populations_->is_gametocytaemic();
}

void GPU::Person::generate_prob_present_at_mda_by_age() {
  if (prob_present_at_mda_by_age().empty()) {
    for (auto i = 0; i < Model::CONFIG->mean_prob_individual_present_at_mda().size(); i++) {
      auto value = Model::RANDOM->random_beta(Model::CONFIG->prob_individual_present_at_mda_distribution()[i].alpha,
                                              Model::CONFIG->prob_individual_present_at_mda_distribution()[i].beta);
      prob_present_at_mda_by_age_.push_back(value);
    }
  }
}

double GPU::Person::prob_present_at_mda() {
  auto i = 0;
  // std::cout << "hello " << i << std::endl;
  while (age_ > Model::CONFIG->age_bracket_prob_individual_present_at_mda()[i]
         && i < Model::CONFIG->age_bracket_prob_individual_present_at_mda().size()) {
    i++;
  }

  return prob_present_at_mda_by_age_[i];
}

bool GPU::Person::has_effective_drug_in_blood() const {
  for (const auto& kv_drug : *drugs_in_blood_->drugs()) {
    if (kv_drug.second->last_update_value() > 0.5) return true;
  }
  return false;
}
double GPU::Person::draw_random_relative_biting_rate(::Random* pRandom, Config* pConfig) {
  auto result =
      pRandom->random_gamma(pConfig->relative_bitting_info().gamma_a, pConfig->relative_bitting_info().gamma_b);

  while (result > (pConfig->relative_bitting_info().max_relative_biting_value
                   - pConfig->relative_bitting_info().min_relative_biting_value)) {
    // re-draw
    result = pRandom->random_gamma(pConfig->relative_bitting_info().gamma_a, pConfig->relative_bitting_info().gamma_b);
  }

  return result + pConfig->relative_bitting_info().min_relative_biting_value;
}

void GPU::Person::generate_render_entity(int location, bool is_circulate){
  float width = Model::CONFIG->debug_config().width > 0 ? Model::CONFIG->debug_config().width : Model::CONFIG->render_config().window_width;
  float height = Model::CONFIG->debug_config().height > 0 ? Model::CONFIG->debug_config().height : Model::CONFIG->render_config().window_height;
  location_col_ = thrust::get<1>(Model::CONFIG->location_db()[location].asc_cell_data);
  location_row_ = thrust::get<2>(Model::CONFIG->location_db()[location].asc_cell_data);
  float unit_x = width/(float)Model::CONFIG->asc_pop_ncols();
  float unit_y = height/(float)Model::CONFIG->asc_pop_nrows();
  float base_x_left = unit_x*location_col_;
  float base_x_right = unit_x*location_col_ + unit_x;
  float base_y_bottom = unit_y*location_row_;
  float base_y_top = unit_y*location_row_+ unit_y;
  float range_x = base_x_right - base_x_left;
  float range_y = base_y_top - base_y_bottom;
  float rand_x = Model::RANDOM->random_uniform_double(0.0,1.0);
  float rand_y = Model::RANDOM->random_uniform_double(0.0,1.0);
  float x = rand_x*range_x + base_x_left;
  float y = height - (rand_y*range_y + base_y_bottom);//OGL from bottom to ptop, so invert Y axis only
  model_ = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, 0.0f));
  color_ = Model::CONFIG->h_location_colors[location];
  if(is_circulate) color_ = glm::vec4(1.0f,1.0f,1.0f,1.0f);
}

void ::GPU::Person::set_liver_parasite_type(GPU::Genotype* value) {
  liver_parasite_type_ = value;
  if(value ==nullptr){
    person_index_gpu->h_person_update_info()[index_].person_liver_parasite_genotype[0] = '\0';
  }
  else{
    std::copy( value->aa_sequence.begin(),
               value->aa_sequence.end(),
               person_index_gpu->h_person_update_info()[index_].person_liver_parasite_genotype);
    person_index_gpu->h_person_update_info()[index_].person_liver_parasite_genotype[MAX_GENOTYPE_LOCUS] = '\0';
  }
}

GPU::Genotype* GPU::Person::liver_parasite_type() const {
  return liver_parasite_type_;
}
