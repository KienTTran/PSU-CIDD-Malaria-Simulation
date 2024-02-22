/* 
 * File:   ImportationEvent.cu
 * Author: Merlin
 * 
 * Created on March 12, 2015, 12:23 PM
 */

#include "ImportationEvent.cuh"
#include "Model.h"
#include "Gpu/Population/Population.cuh"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/Person.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"

GPU::ImportationEvent::
ImportationEvent(const int &location, const int &execute_at, const int &genotype_id, const int &number_of_cases)
    : location_(location),
      genotype_id_(
          genotype_id),
      number_of_cases_(
          number_of_cases) {
  time = execute_at;
}

GPU::ImportationEvent::~ImportationEvent() = default;

void GPU::ImportationEvent::schedule_event(GPU::Scheduler* scheduler, const int &location, const int &execute_at,
                                      const int &genotype_id,
                                      const int &number_of_cases) {
  if (scheduler!=nullptr) {
    auto* e = new GPU::ImportationEvent(location, execute_at, genotype_id, number_of_cases);
    e->dispatcher = nullptr;
    e->time = execute_at;
    scheduler->schedule_population_event(e);

  }
}

void GPU::ImportationEvent::execute() {
  const auto number_of_importation_cases = Model::RANDOM->random_poisson(number_of_cases_);
  auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();

  for (auto i = 0; i < number_of_importation_cases; i++) {

    const size_t ind_ac = Model::RANDOM->random_uniform(static_cast<unsigned long>(pi->vPerson()[0][0].size()));
    if (pi->vPerson()[0][0][ind_ac].empty()) {
      continue;
    }
    const size_t index = Model::RANDOM->random_uniform(pi->vPerson()[0][0][ind_ac].size());
    auto* p = pi->vPerson()[0][0][ind_ac][index];

    p->immune_system()->set_increase(true);
    p->set_host_state(GPU::Person::ASYMPTOMATIC);

    auto* blood_parasite = p->add_new_parasite_to_blood(Model::CONFIG->gpu_genotype_db.at(static_cast<const unsigned long &>(genotype_id_)));

    auto size = Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic;

    blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(size);
    blood_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
    blood_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());

    //        Model::GPU_POPULATION->initial_infection(pi->vPerson()[0][0][ind_ac][index], Model::CONFIG->parasite_db()->get(0));
  }

}
