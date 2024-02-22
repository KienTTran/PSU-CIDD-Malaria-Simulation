/*
 * File:   Population.cu
 * Author: nguyentran
 *
 * Created on April 15, 2013, 10:49 AM
 */

#include "Population.cuh"

#include <cfloat>
#include <cmath>

#include "Constants.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Events/BirthdayEvent.cuh"
#include "Gpu/Events/SwitchImmuneComponentEvent.cuh"
#include "Helpers/ObjectHelpers.h"
#include "Helpers/TimeHelpers.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Model.h"
#include "Properties/PersonIndexAll.cuh"
#include "Spatial/SpatialModel.hxx"
#include "Helpers/UniqueId.hxx"
#include "easylogging++.h"
#include "Gpu/Mosquito/Mosquito.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationMovingLevel.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/InfantImmuneComponent.cuh"
#include "Gpu/Population/NonInfantImmuneComponent.cuh"
#include "Gpu/Population/Person.cuh"
#include "Gpu/Population/PopulationKernel.cuh"
#include "Gpu/Core/Scheduler.cuh"


GPU::Population::Population(Model* model) : model_(model) {
  person_index_list_ = new GPUPersonIndexPtrList();
  all_persons_ = new GPU::PersonIndexAll();

  person_index_list_->push_back(all_persons_);
}
GPU::Population::~Population() {
  // release memory for all persons
  if (all_persons_ != nullptr) {
    for (auto& person : all_persons_->vPerson()) {
      ObjectHelpers::delete_pointer<GPU::Person>(person);
    }
    all_persons_->vPerson().clear();
    all_persons_ = nullptr;
  }

  // release person_indexes

  if (person_index_list_ != nullptr) {
    for (GPU::PersonIndex* person_index : *person_index_list_) {
      ObjectHelpers::delete_pointer<GPU::PersonIndex>(person_index);
    }

    person_index_list_->clear();
    ObjectHelpers::delete_pointer<GPUPersonIndexPtrList>(person_index_list_);
  }
}

void GPU::Population::add_person(GPU::Person* person) {
  for (GPU::PersonIndex* person_index : *person_index_list_) {
    person_index->add(person);
  }
  person->set_population(this);
}

void GPU::Population::remove_person(GPU::Person* person) {
  for (GPU::PersonIndex* person_index : *person_index_list_) {
    person_index->remove(person);
  }
  person->set_population(nullptr);
}

void GPU::Population::remove_dead_person(GPU::Person* person) {
  remove_person(person);
  ObjectHelpers::delete_pointer<GPU::Person>(person);
}

void GPU::Population::notify_change(GPU::Person* p, const GPU::Person::Property& property, const void* oldValue,
                               const void* newValue) {
  for (GPU::PersonIndex* person_index : *person_index_list_) {
    person_index->notify_change(p, property, oldValue, newValue);
  }
}

std::size_t GPU::Population::size(const int& location, const int& age_class) {
  if (location == -1) {
    return all_persons_->size();
  }
  auto pi_lsa = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();

  if (pi_lsa == nullptr) {
    return 0;
  }
  std::size_t temp = 0;
  if (age_class == -1) {
    for (auto state = 0; state < GPU::Person::NUMBER_OF_STATE - 1; state++) {
      for (auto ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        temp += pi_lsa->vPerson()[location][state][ac].size();
      }
    }
  } else {
    for (auto state = 0; state < GPU::Person::NUMBER_OF_STATE - 1; state++) {
      temp += pi_lsa->vPerson()[location][state][age_class].size();
    }
  }
  return temp;
}

std::size_t GPU::Population::size(const int& location, const GPU::Person::HostStates& hs, const int& age_class) {
  if (location == -1) {
    return all_persons_->size();
  }
  auto pi_lsa = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  return (pi_lsa->vPerson()[location][hs][age_class].size());
}

// new
std::size_t GPU::Population::size_residents_only(const int& location) {
  if (location == -1) {
    return all_persons_->size();
  }

  auto pi_lsa = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();

  if (pi_lsa == nullptr) {
    return 0;
  }
  auto temp = 0ul;
  for (auto state = 0; state < GPU::Person::NUMBER_OF_STATE - 1; state++) {
    for (auto ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
      for (auto i = 0; i < pi_lsa->vPerson()[location][state][ac].size(); i++) {
        if (pi_lsa->vPerson()[location][state][ac][i]->residence_location() == location) {
          temp++;
        }
      }
    }
  }
  return temp;
}

void GPU::Population::perform_infection_event() {
  auto start = std::chrono::high_resolution_clock::now();
  //    std::cout << "Infection Event" << std::endl;

  GPUPersonPtrVector today_infections;
  auto tracking_index = Model::GPU_SCHEDULER->current_time() % Model::CONFIG->number_of_tracking_days();
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    const auto force_of_infection = force_of_infection_for_N_days_by_location[tracking_index][loc];
    if (force_of_infection <= DBL_EPSILON) continue;

    const auto new_beta = Model::CONFIG->location_db()[loc].beta
                          * Model::CONFIG->seasonal_info()->get_seasonal_factor(Model::GPU_SCHEDULER->calendar_date, loc);

    auto poisson_means = new_beta * force_of_infection;

    auto number_of_bites = Model::RANDOM->random_poisson(poisson_means);
    if (number_of_bites <= 0) continue;

    // data_collector store number of bites
    Model::GPU_DATA_COLLECTOR->collect_number_of_bites(loc, number_of_bites);

    auto persons_bitten_today = Model::RANDOM->roulette_sampling<GPU::Person>(
        number_of_bites, individual_relative_biting_by_location[loc], all_alive_persons_by_location[loc], false,
        sum_relative_biting_by_location[loc]);

    for (auto* person : persons_bitten_today) {
      assert(person->host_state() != GPU::Person::DEAD);
      person->increase_number_of_times_bitten();

      auto genotype_id = Model::GPU_MOSQUITO->random_genotype(loc, tracking_index);

      if(genotype_id == -1) continue;

      const auto p_infectious = Model::RANDOM->random_flat(0.0, 1.0);
      // only infect with real infectious bite
      if (Model::CONFIG->using_variable_probability_infectious_bites_cause_infection()) {
        if (p_infectious <= person->p_infection_from_an_infectious_bite()) {
          if (person->host_state() != GPU::Person::EXPOSED && person->liver_parasite_type() == nullptr) {
            person->today_infections()->push_back(genotype_id);
            today_infections.push_back(person);
          }
        }
      } else if (p_infectious <= Model::CONFIG->p_infection_from_an_infectious_bite()) {
        if (person->host_state() != GPU::Person::EXPOSED && person->liver_parasite_type() == nullptr) {
          person->today_infections()->push_back(genotype_id);
          today_infections.push_back(person);
        }
      }
    }
  }
  //    std::cout << "Solve infections"<< std::endl;
  // solve Multiple infections
  if (today_infections.empty()) return;

  for (auto* p : today_infections) {
    if (!p->today_infections()->empty()) {
      Model::GPU_DATA_COLLECTOR->monthly_number_of_new_infections_by_location()[p->location()] += 1;
    }
    p->randomly_choose_parasite();
  }

  today_infections.clear();
  auto lapse = std::chrono::high_resolution_clock::now() - start;
  if(Model::CONFIG->debug_config().enable_debug_text){
      LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[Population] Update population infection CPU time: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
}



void GPU::Population::initialize() {
  if (model() != nullptr) {
    // those vector will be used in the initial infection
    const auto number_of_locations = Model::CONFIG->number_of_locations();

    individual_relative_biting_by_location =
        std::vector<std::vector<double>>(number_of_locations, std::vector<double>());
    individual_relative_moving_by_location =
        std::vector<std::vector<double>>(number_of_locations, std::vector<double>());
    individual_foi_by_location = std::vector<std::vector<double>>(number_of_locations, std::vector<double>());

    all_alive_persons_by_location = std::vector<std::vector<GPU::Person*>>(number_of_locations, std::vector<GPU::Person*>());

    sum_relative_biting_by_location = std::vector<double>(number_of_locations, 0);
    sum_relative_moving_by_location = std::vector<double>(number_of_locations, 0);

    current_force_of_infection_by_location = std::vector<double>(number_of_locations, 0);

    force_of_infection_for_N_days_by_location = std::vector<std::vector<double>>(
        Model::CONFIG->number_of_tracking_days(), std::vector<double>(number_of_locations, 0));

    // initalize other person index
    initialize_person_indices();

    // initialize population
    for (auto loc = 0; loc < number_of_locations; loc++) {
      const auto popsize_by_location = static_cast<int>(Model::CONFIG->location_db()[loc].population_size
                                                        * Model::CONFIG->artificial_rescaling_of_population_size());
      auto temp_sum = 0;
      for (auto age_class = 0; age_class < Model::CONFIG->initial_age_structure().size(); age_class++) {
        int number_of_individual_by_loc_age_class = static_cast<int>(popsize_by_location
                                                                     * Model::CONFIG->location_db()[loc].age_distribution[age_class]);
//            printf("[Population] Init loc %d ac %d (%d + %d)",loc, age_class, number_of_individual_by_loc_age_class,temp_sum);
        if(Model::CONFIG->location_db()[loc].population_size > temp_sum + number_of_individual_by_loc_age_class){
          temp_sum += number_of_individual_by_loc_age_class;
          if(age_class == Model::CONFIG->initial_age_structure().size() - 1) {
            number_of_individual_by_loc_age_class = Model::CONFIG->location_db()[loc].population_size - temp_sum;
          }
        }
        else{
          if(age_class == Model::CONFIG->initial_age_structure().size() - 1){
            number_of_individual_by_loc_age_class = Model::CONFIG->location_db()[loc].population_size - temp_sum;
            if(number_of_individual_by_loc_age_class < 0){
              number_of_individual_by_loc_age_class = 0;
            }
          }
          else{
            int minus = (temp_sum + number_of_individual_by_loc_age_class) - Model::CONFIG->location_db()[loc].population_size;
            number_of_individual_by_loc_age_class = number_of_individual_by_loc_age_class - minus;
            temp_sum += number_of_individual_by_loc_age_class;
          }
        }
        for (auto i = 0; i < number_of_individual_by_loc_age_class; i++) {
            generate_individual(loc, age_class);
        }
      }
    }
  }
}

void GPU::Population::introduce_initial_cases() {
  if (model_ != nullptr) {
    auto start = std::chrono::high_resolution_clock::now();
    // std::cout << Model::CONFIG->initial_parasite_info().size() << std::endl;
    for (const auto p_info : Model::CONFIG->initial_parasite_info()) {
      auto num_of_infections = Model::RANDOM->random_poisson(std::round(size(p_info.location) * p_info.prevalence));
      num_of_infections = num_of_infections <= 0 ? 1 : num_of_infections;

      auto* genotype = Model::CONFIG->gpu_genotype_db.at(p_info.parasite_type_id);
//       std::cout << p_info.location << "-" << p_info.parasite_type_id << "-" << num_of_infections << std::endl;
      introduce_parasite(p_info.location, genotype, num_of_infections);
    }
    // update current foi
    update_current_foi();
//    Model::GPU_POPULATION_KERNEL->update_current_foi();

    // update force of infection for N days
    for (auto d = 0; d < Model::CONFIG->number_of_tracking_days(); d++) {
      for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        force_of_infection_for_N_days_by_location[d][loc] = current_force_of_infection_by_location[loc];
      }
      Model::GPU_MOSQUITO->infect_new_cohort_in_PRMC(Model::CONFIG, Model::RANDOM, this, d);
    }
    auto lapse = std::chrono::high_resolution_clock::now() - start;
    LOG(INFO) << fmt::format("[Population] introduce_parasite CPU time: {} ms",
                             std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
  }
}

void GPU::Population::introduce_parasite(const int& location, GPU::Genotype* parasite_type, const int& num_of_infections) {
  if (model_ != nullptr) {
    auto persons_bitten_today = Model::RANDOM->roulette_sampling<GPU::Person>(
        num_of_infections, Model::GPU_POPULATION->individual_relative_biting_by_location[location],
        Model::GPU_POPULATION->all_alive_persons_by_location[location], false);

    for (auto* person : persons_bitten_today) {
      initial_infection(person, parasite_type);
    }
  }
}

void GPU::Population::initial_infection(GPU::Person* person, GPU::Genotype* parasite_type) const {
  if(person == nullptr) return;

  auto* blood_parasite = person->add_new_parasite_to_blood(parasite_type);

  const auto size =
      model_->random()->random_flat(Model::CONFIG->parasite_density_level().log_parasite_density_from_liver,
                                    Model::CONFIG->parasite_density_level().log_parasite_density_clinical);

  blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
  blood_parasite->set_last_update_log10_parasite_density(size);

  const auto p_clinical = person->get_probability_progress_to_clinical();
  const auto p = model_->random()->random_flat(0.0, 1.0);

  if (p < p_clinical) {
    // progress to clinical after several days
    blood_parasite->set_update_function(model_->gpu_progress_to_clinical_update_function());
    person->schedule_progress_to_clinical_event_by(blood_parasite);
  } else {
    // only progress to clearance by Immune system
    // progress to clearance
    blood_parasite->set_update_function(model_->gpu_immunity_clearance_update_function());
  }

//  assert(person->id() == pi->h_persons()[p_index]->id());
//  assert(person->index() == pi->h_persons()[p_index]->index());

    for(auto *parasite: *person->all_clonal_parasite_populations()->parasites()){
        LOG_IF(person->index() >= 1040 && person->index() <= 1045,INFO)
            << fmt::format("{} GPU::Population::initial_infection {} {} {} {}",
               person->index(),
               parasite->index(),
               parasite->update_function()->type(),
               parasite->genotype()->aa_sequence.c_str(),
               parasite->last_update_log10_parasite_density());
    }
}

void GPU::Population::perform_birth_event() {
  //    std::cout << "Birth Event" << std::endl;
  auto tp_start = std::chrono::high_resolution_clock::now();
  int birth_sum = 0;

  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    auto poisson_means = size(loc) * Model::CONFIG->birth_rate() / Constants::DAYS_IN_YEAR();
    const auto number_of_births = Model::RANDOM->random_poisson(poisson_means);
    for (auto i = 0; i < number_of_births; i++) {
      give_1_birth(loc);
      Model::GPU_DATA_COLLECTOR->update_person_days_by_years(loc, Constants::DAYS_IN_YEAR() - Model::GPU_SCHEDULER->current_day_in_year());
      birth_sum++;
    }
    Model::CONFIG->location_db()[loc].population_size += number_of_births;
  }

  if(Model::CONFIG->debug_config().enable_debug_text){
    auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
    auto *pi = get_person_index<GPU::PersonIndexGPU>();
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
    << "[Population] Update population birth (" << birth_sum << " " << pi->h_persons().size() << ") CPU time: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
}


void GPU::Population::generate_individual(int location, int age_class) {
  auto p = new GPU::Person();
  p->init();
  p->set_id(UniqueId::get_instance().get_uid());
  p->set_location(location);
  p->set_residence_location(location);
  p->set_host_state(GPU::Person::SUSCEPTIBLE);

  const auto age_from = (age_class == 0) ? 0 : Model::CONFIG->initial_age_structure()[age_class - 1];
  const auto age_to = Model::CONFIG->initial_age_structure()[age_class];

  // std::cout << i << "\t" << age_class << "\t" << age_from << "\t" << age_to << std::endl;

  // set age will also set ageclass
  p->set_age(static_cast<const int&>(Model::RANDOM->random_uniform_int(age_from, age_to + 1)));
  //                    std::cout << p->age() << " \t" << p->age_class() << std::endl;
  //                    p->set_age_class(age_class);

  int days_to_next_birthday = Model::RANDOM->random_uniform(Constants::DAYS_IN_YEAR());

  auto simulation_time_birthday = TimeHelpers::get_simulation_time_birthday(days_to_next_birthday, p->age(),
                                                                            Model::GPU_SCHEDULER->calendar_date);
  p->set_birthday(simulation_time_birthday);

  LOG_IF(simulation_time_birthday > 0, FATAL)
    << "simulation_time_birthday have to be <= 0 when initilizing population";

  GPU::BirthdayEvent::schedule_event(Model::GPU_SCHEDULER, p, days_to_next_birthday);

  // set immune component
  if (simulation_time_birthday + Constants::DAYS_IN_YEAR() / 2 >= 0) {
    LOG_IF(p->age() > 0, FATAL) << "Error in calculating simulation_time_birthday";
    // LOG(INFO) << "Infant: " << p->age() << " - " << simulation_time_birthday;
    p->immune_system()->set_immune_component(new GPU::InfantImmuneComponent());
    // schedule for switch
    GPU::SwitchImmuneComponentEvent::schedule_for_switch_immune_component_event(
            Model::GPU_SCHEDULER, p, simulation_time_birthday + Constants::DAYS_IN_YEAR() / 2);
  } else {
    // LOG(INFO) << "Adult: " << p->age() << " - " << simulation_time_birthday;
    p->immune_system()->set_immune_component(new GPU::NonInfantImmuneComponent());
  }

  auto immune_value = Model::RANDOM->random_beta(Model::CONFIG->immune_system_information().alpha_immune,
                                                 Model::CONFIG->immune_system_information().beta_immune);
  p->immune_system()->set_latest_immune_value(immune_value);
  p->immune_system()->set_increase(false);
  //                    p->draw_random_immune();

  p->innate_relative_biting_rate = GPU::Person::draw_random_relative_biting_rate(Model::RANDOM, Model::CONFIG);
  p->update_relative_bitting_rate();

  p->set_moving_level(Model::CONFIG->moving_level_generator().draw_random_level(Model::RANDOM));

  p->set_latest_update_time(0);

  int time = Model::RANDOM->random_uniform(Model::CONFIG->update_frequency()) + 1;
  p->generate_prob_present_at_mda_by_age();

  add_person(p);

  individual_relative_biting_by_location[location].push_back(p->current_relative_biting_rate);
  individual_relative_moving_by_location[location].push_back(
          Model::CONFIG->circulation_info().v_moving_level_value[p->moving_level()]);

  sum_relative_biting_by_location[location] += p->current_relative_biting_rate;
  sum_relative_moving_by_location[location] +=
          Model::CONFIG->circulation_info().v_moving_level_value[p->moving_level()];

  all_alive_persons_by_location[location].push_back(p);
}

void GPU::Population::give_1_birth(const int& location) {
  auto p = new GPU::Person();
  p->init();
  p->set_id(UniqueId::get_instance().get_uid());
  p->set_age(0);
  p->set_host_state(GPU::Person::SUSCEPTIBLE);
  p->set_age_class(0);
  p->set_location(location);
  p->set_residence_location(location);
  p->immune_system()->set_immune_component(new GPU::InfantImmuneComponent());
  p->immune_system()->set_latest_immune_value(1.0);
  p->immune_system()->set_increase(false);

  p->set_latest_update_time(Model::GPU_SCHEDULER->current_time());

//  p->draw_random_immune();

  // set_relative_biting_rate
  p->innate_relative_biting_rate = GPU::Person::draw_random_relative_biting_rate(Model::RANDOM, Model::CONFIG);
  p->update_relative_bitting_rate();

  p->set_moving_level(Model::CONFIG->moving_level_generator().draw_random_level(Model::RANDOM));

  p->set_birthday(Model::GPU_SCHEDULER->current_time());
  const auto number_of_days_to_next_birthday =
      TimeHelpers::number_of_days_to_next_year(Model::GPU_SCHEDULER->calendar_date);
  GPU::BirthdayEvent::schedule_event(Model::GPU_SCHEDULER, p,
                                Model::GPU_SCHEDULER->current_time() + number_of_days_to_next_birthday);

  // schedule for switch
  GPU::SwitchImmuneComponentEvent::schedule_for_switch_immune_component_event(
      Model::GPU_SCHEDULER, p, Model::GPU_SCHEDULER->current_time() + Constants::DAYS_IN_YEAR() / 2);

  //    p->startLivingTime = (Global::startTreatmentDay > Global::scheduler->currentTime) ? Global::startTreatmentDay :
  //    Global::scheduler->currentTime;
  p->generate_prob_present_at_mda_by_age();

  add_person(p);
}

void GPU::Population::perform_death_event() {
  auto tp_start = std::chrono::high_resolution_clock::now();
//      std::cout << "Death Event" << std::endl;
  // simply change state to dead and release later
  auto pi = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  if (pi == nullptr) return;
  auto& location_db = Model::CONFIG->location_db();

  int dead_sum = 0;
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    int loc_deaths = 0;
    for (auto hs = 0; hs < GPU::Person::NUMBER_OF_STATE - 1; hs++) {
      if (hs == GPU::Person::DEAD) continue;
      for (auto ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        const size_t size = pi->vPerson()[loc][hs][ac].size();
        if (size == 0) continue;
        auto poisson_means = size * Model::CONFIG->death_rate_by_age_class()[ac] / Constants::DAYS_IN_YEAR();

        assert(Model::CONFIG->death_rate_by_age_class().size() == Model::CONFIG->number_of_age_classes());
        const auto number_of_deaths = Model::RANDOM->random_poisson(poisson_means);
        if (number_of_deaths == 0) continue;
        for (int i = 0; i < number_of_deaths; i++) {
          // change state to Death;
          int index = Model::RANDOM->random_uniform(pi->vPerson()[loc][hs][ac].size());
//          LOG_IF(loc == 0 && hs == 2 && ac == 14,INFO)
//                  << fmt::format("[Population] Death before select loc {} hs {} ac {} deaths {} size {}-{} size2 {} size 3 {} index {}",
//                                 loc, hs, ac, number_of_deaths, size, pi->vPerson()[loc][hs][ac].size(),
//                                 pi2->h_person_host_states().size(), pi2->h_persons().size(), index);
          auto* p = pi->vPerson()[loc][hs][ac][index];
          p->cancel_all_events_except(nullptr);
//          LOG_IF(loc == 0 && hs == 2 && ac == 14,INFO)
//                  << fmt::format("[Population] Death before set host state loc {} hs {} ac {} deaths {} size {}-{} size2 {} size 3 {} index {}",
//                                 loc, hs, ac, number_of_deaths, size,pi->vPerson()[loc][hs][ac].size(),
//                                 pi2->h_person_host_states().size(), pi2->h_persons().size(), index);
          p->set_host_state(GPU::Person::DEAD);
//          LOG_IF(loc == 0 && hs == 2 && ac == 14,INFO)
//                  << fmt::format("[Population] Death after set host state loc {} hs {} ac {} deaths {} size {}-{} size2 {} size 3 {} index {}",
//                                 loc, hs, ac, number_of_deaths, size, pi->vPerson()[loc][hs][ac].size(),
//                                 pi2->h_person_host_states().size(), pi2->h_persons().size(),index);
        }
        loc_deaths += number_of_deaths;
        dead_sum += number_of_deaths;
      }
    }
    location_db[loc].population_size -= loc_deaths;
  }
  clear_all_dead_state_individual();
  if(Model::CONFIG->debug_config().enable_debug_text){
    auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
    auto *pi = get_person_index<GPU::PersonIndexGPU>();
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
    << "[Population] Update population death (" << dead_sum << " " << pi->h_persons().size() << ") CPU time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
}

void GPU::Population::clear_all_dead_state_individual() {
  // return all Death to object pool and clear vPersonIndex[l][dead][ac] for all location and ac
  auto pi = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  GPUPersonPtrVector removePersons;

  for (int loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
      for (auto person : pi->vPerson()[loc][GPU::Person::DEAD][ac]) {
        removePersons.push_back(person);
      }
    }
  }

  for (GPU::Person* p : removePersons) {
    remove_dead_person(p);
  }
}

void GPU::Population::perform_circulation_event() {
  auto tp_start = std::chrono::high_resolution_clock::now();
  // for each location
  //  get number of circulations based on size * circulation_percent
  //  distributes that number into others location based of other location size
  //  for each number in that list select an individual, and schedule a movement event on next day
  GPUPersonPtrVector today_circulations;

  std::vector<int> v_number_of_residents_by_location(Model::CONFIG->number_of_locations(), 0);

  for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    //        v_number_of_residents_by_location[target_location] = (size(target_location));
    v_number_of_residents_by_location[location] = Model::GPU_DATA_COLLECTOR->popsize_residence_by_location()[location];
    //        std::cout << v_original_pop_size_by_location[target_location] << std::endl;
  }

  for (int from_location = 0; from_location < Model::CONFIG->number_of_locations(); from_location++) {
    auto poisson_means = size(from_location) * Model::CONFIG->circulation_info().circulation_percent;
//    LOG_IF(poisson_means == 0, DEBUG)
//      << "[Population] Update population circulation CPU " << from_location << " "
//      << Model::GPU_DATA_COLLECTOR->popsize_residence_by_location()[from_location]  << " poisson_means = 0";
    if (poisson_means == 0) continue;
    const auto number_of_circulating_from_this_location = Model::RANDOM->random_poisson(poisson_means);
//    LOG_IF(number_of_circulating_from_this_location == 0, DEBUG)
//      << "[Population] Update population circulation CPU " << from_location << " "
//      << Model::GPU_DATA_COLLECTOR->popsize_residence_by_location()[from_location]  << " number_of_circulating_from_this_location = 0";
    if (number_of_circulating_from_this_location == 0) continue;

    DoubleVector v_relative_outmovement_to_destination(Model::CONFIG->number_of_locations(), 0);
    v_relative_outmovement_to_destination = Model::CONFIG->spatial_model()->get_v_relative_out_movement_to_destination(
        from_location, Model::CONFIG->number_of_locations(), Model::CONFIG->spatial_distance_matrix()[from_location],
        v_number_of_residents_by_location);

    std::vector<unsigned int> v_num_leavers_to_destination(
        static_cast<unsigned int>(Model::CONFIG->number_of_locations()));

    Model::RANDOM->random_multinomial(static_cast<int>(v_relative_outmovement_to_destination.size()),
                                      static_cast<unsigned int>(number_of_circulating_from_this_location),
                                      &v_relative_outmovement_to_destination[0], &v_num_leavers_to_destination[0]);

    for (int target_location = 0; target_location < Model::CONFIG->number_of_locations(); target_location++) {
      //            std::cout << v_num_leavers_to_destination[target_location] << std::endl;
      if (v_num_leavers_to_destination[target_location] == 0) continue;
      //            std::cout << Model::GPU_SCHEDULER->current_time() << "\t" << from_location << "\t" << target_location <<
      //            "\t"
      //                      << v_num_leavers_to_destination[target_location] << std::endl;
      perform_circulation_for_1_location(from_location, target_location, v_num_leavers_to_destination[target_location],
                                         today_circulations);
    }
  }

  for (auto* p : today_circulations) {
    p->randomly_choose_target_location();
  }

  auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
  if(Model::CONFIG->debug_config().enable_debug_text){
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
    << "[Population] Update population circulation CPU (" << today_circulations.size() << ") event time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
  today_circulations.clear();
}

void GPU::Population::perform_circulation_for_1_location(const int& from_location, const int& target_location,
                                                    const int& number_of_circulations,
                                                    std::vector<GPU::Person*>& today_circulations) {
  auto persons_moving_today = Model::RANDOM->roulette_sampling<GPU::Person>(
      number_of_circulations, individual_relative_moving_by_location[from_location],
      all_alive_persons_by_location[from_location], false, sum_relative_moving_by_location[from_location]);

  for (auto* person : persons_moving_today) {
    assert(person->host_state() != GPU::Person::DEAD);

    person->today_target_locations()->push_back(target_location);
    today_circulations.push_back(person);
  }
}

bool GPU::Population::has_0_case() {
  auto pi = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  for (int loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    for (int hs = GPU::Person::EXPOSED; hs <= GPU::Person::CLINICAL; hs++) {
      for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        if (!pi->vPerson()[loc][hs][ac].empty()) {
          return false;
        }
      }
    }
  }
  return true;
}

void GPU::Population::initialize_person_indices() {
  const int number_of_location = Model::CONFIG->number_of_locations();
  const int number_of_hoststate = GPU::Person::NUMBER_OF_STATE;
  const int number_of_ageclasses = Model::CONFIG->number_of_age_classes();

  auto p_index_by_l_s_a =
      new GPU::PersonIndexByLocationStateAgeClass(number_of_location, number_of_hoststate, number_of_ageclasses);
  person_index_list_->push_back(p_index_by_l_s_a);

  auto p_index_location_moving_level = new GPU::PersonIndexByLocationMovingLevel(
      number_of_location, Model::CONFIG->circulation_info().number_of_moving_levels);
    person_index_list_->push_back(p_index_location_moving_level);

    auto p_index_gpu = new GPU::PersonIndexGPU();
    person_index_list_->push_back(p_index_gpu);
}

void GPU::Population::update_all_individuals() {
  auto start = std::chrono::high_resolution_clock::now();
  // update all individuals
  auto pi = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  for (int loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    for (int hs = 0; hs < Person::DEAD; hs++) {
      for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        for (auto* person : pi->vPerson()[loc][hs][ac]) {
          person->update();
        }
      }
    }
  }
  auto lapse = std::chrono::high_resolution_clock::now() - start;
  if(Model::CONFIG->debug_config().enable_debug_text){
      LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[Population] Update population all individuals CPU time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
//  auto *pi2 = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
//  for(int i = pi2->h_persons().size() - 5; i < pi2->h_persons().size(); i++){
//    LOG_IF(Model::CONFIG->debug_config().enable_debug_text,INFO)
//      << fmt::format("[PopulationKernel update_all_individuals] CPU Person {} last update time {} parasites_size {}"
//                     " person_biting {} person_moving {}",
//                     i,
//                     pi2->h_persons()[i]->latest_update_time(),
//                     pi2->h_persons()[i]->all_clonal_parasite_populations()->size(),
//                     pi2->h_persons()[i]->current_relative_biting_rate,
//                     Model::CONFIG->circulation_info().v_moving_level_value[pi2->h_persons()[i]->moving_level()]);
//  }
}

void GPU::Population::persist_current_force_of_infection_to_use_N_days_later() {
  auto start = std::chrono::high_resolution_clock::now();
  for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    force_of_infection_for_N_days_by_location[Model::GPU_SCHEDULER->current_time()
                                              % Model::CONFIG->number_of_tracking_days()][loc] =
        current_force_of_infection_by_location[loc];
  }
  auto lapse = std::chrono::high_resolution_clock::now() - start;
  if(Model::CONFIG->debug_config().enable_debug_text){
    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[Population] Update FOI N days CPU time: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
}

void GPU::Population::update_current_foi() {
  auto start = std::chrono::high_resolution_clock::now();
  auto pi = get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  for (int loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
    // reset force of infection for each location
    current_force_of_infection_by_location[loc] = 0.0;
    sum_relative_biting_by_location[loc] = 0.0;
    sum_relative_moving_by_location[loc] = 0.0;

    // using clear so as system will not reallocate memory slot for vector
    individual_foi_by_location[loc].clear();
    individual_relative_biting_by_location[loc].clear();
    individual_relative_moving_by_location[loc].clear();
    all_alive_persons_by_location[loc].clear();

    for (int hs = 0; hs < GPU::Person::DEAD; hs++) {
      for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        for (auto* person : pi->vPerson()[loc][hs][ac]) {
          double log_10_total_infectious_density =
              person->all_clonal_parasite_populations()->log10_total_infectious_denstiy;

          auto individual_foi = log_10_total_infectious_density == GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY
                                    ? 0.0
                                    : person->current_relative_biting_rate
                                          * GPU::Person::relative_infectivity(log_10_total_infectious_density);
          individual_foi_by_location[loc].push_back(individual_foi);
          individual_relative_biting_by_location[loc].push_back(person->current_relative_biting_rate);
          individual_relative_moving_by_location[loc].push_back(
              Model::CONFIG->circulation_info().v_moving_level_value[person->moving_level()]);

          sum_relative_biting_by_location[loc] += person->current_relative_biting_rate;
          sum_relative_moving_by_location[loc] += Model::CONFIG->circulation_info().v_moving_level_value[person->moving_level()];
          current_force_of_infection_by_location[loc] += individual_foi;
          all_alive_persons_by_location[loc].push_back(person);
        }
      }
    }
  }
  auto lapse = std::chrono::high_resolution_clock::now() - start;
  if(Model::CONFIG->debug_config().enable_debug_text){
      LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
      << "[Population] Update population current foi CPU time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count() << " ms";
  }
//  if(Model::CONFIG->debug_config().enable_debug_text){
//    for (int loc = Model::CONFIG->number_of_locations() - 5; loc < Model::CONFIG->number_of_locations(); loc++) {
//      printf("%d CPU sum_relative_biting_by_location[%d] biting %f\n",Model::GPU_SCHEDULER->current_time(),loc,sum_relative_biting_by_location[loc]);
//      printf("%d CPU sum_relative_moving_by_location[%d] moving %f\n",Model::GPU_SCHEDULER->current_time(),loc,sum_relative_moving_by_location[loc]);
//      printf("%d CPU sum_relative_moving_by_location[%d] foi %f\n",Model::GPU_SCHEDULER->current_time(),loc,current_force_of_infection_by_location[loc]);
//    }
//  }
}