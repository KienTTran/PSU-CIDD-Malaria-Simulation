/*
 * File:   ImportationEvent.cu
 * Author: Merlin
 *
 * Created on February 21, 2014, 2:42 PM
 */

#include "ImportationPeriodicallyEvent.cuh"

#include <easylogging++.h>

#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Model.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"

GPU::ImportationPeriodicallyEvent::ImportationPeriodicallyEvent(const int& location, const int& duration, int genotype_id,
                                                           const int& number_of_cases, const int& start_day)
    : location_(location), duration_(duration), genotype_id_(genotype_id), number_of_cases_(number_of_cases) {
  // TODO: remove start_day_
  time = start_day;
}

GPU::ImportationPeriodicallyEvent::~ImportationPeriodicallyEvent() = default;

void GPU::ImportationPeriodicallyEvent::schedule_event(GPU::Scheduler* scheduler, const int& location, const int& duration,
                                                  unsigned int genotype_id, const int& number_of_cases,
                                                  const int& start_day) {
  if (scheduler != nullptr) {
    auto* e = new GPU::ImportationPeriodicallyEvent(location, duration, genotype_id, number_of_cases, start_day);
    e->dispatcher = nullptr;
    e->time = start_day;
    scheduler->schedule_population_event(e);
  }
}

void GPU::ImportationPeriodicallyEvent::execute() {
  // std::cout << date::year_month_day{ Model::GPU_SCHEDULER->calendar_date } << ":import periodically event" << std::endl;
  // schedule importation for the next day
  schedule_event(Model::GPU_SCHEDULER, location_, duration_, genotype_id_, number_of_cases_,
                 Model::GPU_SCHEDULER->current_time() + 1);

  const auto number_of_importation_cases =
      Model::RANDOM->random_poisson(static_cast<double>(number_of_cases_) / duration_);
  if (Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[location_][0] < number_of_importation_cases) {
    return;
  }

  //    std::cout << number_of_cases_ << std::endl;
  auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();
  VLOG_IF(number_of_importation_cases > 0, 2)
      << "Day: " << Model::GPU_SCHEDULER->current_time() << " - Importing " << number_of_importation_cases
      << " at location " << location_ << " with genotype " << genotype_id_;
  for (auto i = 0; i < number_of_importation_cases; i++) {
    std::size_t ind_ac = Model::RANDOM->random_uniform(static_cast<unsigned long>(pi->vPerson()[location_][0].size()));
    if (pi->vPerson()[location_][0][ind_ac].empty()) {
      continue;
    }

    std::size_t index = Model::RANDOM->random_uniform(pi->vPerson()[location_][0][ind_ac].size());
    auto* p = pi->vPerson()[location_][0][ind_ac][index];

    p->immune_system()->set_increase(true);
    p->set_host_state(GPU::Person::ASYMPTOMATIC);

    // check and draw random Genotype
    GPU::Genotype* imported_genotype = nullptr;

    // TODO: rework on this function to have a random genotype string
    uint32_t random_id = Model::RANDOM->random_uniform_int(0, 1);

    switch (genotype_id_) {
      case -1:
        //      new genotype will have 50% chance of 580Y and 50% plasmepsin-2 copy, last allele will always be x
        if (random_id % 2 == 1) {
          random_id -= 1;
        }
        imported_genotype = Model::CONFIG->gpu_genotype_db.at(random_id);
        break;
      case -2:
        // all random even last xX locus new genotype will have
        // 50% chance of 580Y and 50% plasmepsin-2 copy and %50 X ....
        imported_genotype = Model::CONFIG->gpu_genotype_db.at(random_id);
        break;
      default:
        imported_genotype = Model::CONFIG->gpu_genotype_db.at(genotype_id_);
    }

    auto* blood_parasite = p->add_new_parasite_to_blood(imported_genotype);
    //    std::cout << "hello"<< std::endl;

    auto size = Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic;

    blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(size);
    blood_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());

    //        Model::GPU_POPULATION->initial_infection(pi->vPerson()[0][0][ind_ac][index],
    //        Model::CONFIG->parasite_db()->get(0));
  }
}