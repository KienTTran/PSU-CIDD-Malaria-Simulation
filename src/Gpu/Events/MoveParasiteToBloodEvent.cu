/* 
 * File:   MoveParasiteToBloodEvent.cu
 * Author: Merlin
 * 
 * Created on July 31, 2013, 11:14 PM
 */

#include "MoveParasiteToBloodEvent.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Parasites/Genotype.cuh"
#include "Gpu/Population/DrugsInBlood.cuh"
#include "Gpu/Population/Person.cuh"

GPU::MoveParasiteToBloodEvent::MoveParasiteToBloodEvent() : infection_genotype_(nullptr) {}

GPU::MoveParasiteToBloodEvent::~MoveParasiteToBloodEvent() {}

void
GPU::MoveParasiteToBloodEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::Genotype *infection_type, const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new GPU::MoveParasiteToBloodEvent();
    e->dispatcher = p;
    e->set_infection_genotype(infection_type);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::MoveParasiteToBloodEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  LOG_IF(person->index() >= 1040 && person->index() <= 1045,INFO)
    << fmt::format("GPU::MoveParasiteToBloodEvent::execute() {}",person->index());
  auto *parasite_type = person->liver_parasite_type();
  person->set_liver_parasite_type(nullptr);

  //add to blood
  if (person->host_state()==GPU::Person::EXPOSED) {
    person->set_host_state(GPU::Person::ASYMPTOMATIC);
  }

  person->immune_system()->set_increase(true);

  auto new_parasite = person->add_new_parasite_to_blood(parasite_type);

  new_parasite->set_last_update_log10_parasite_density(
      Model::RANDOM->random_normal_truncated(
          Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic, 0.5));

  if (person->has_effective_drug_in_blood()) {
    //person has drug in blood
    new_parasite->set_update_function(Model::MODEL->gpu_having_drug_update_function());
  } else {

    if (person->all_clonal_parasite_populations()->size() > 1) {
      if (Model::CONFIG->allow_new_coinfection_to_cause_symtoms()) {
        person->determine_clinical_or_not(new_parasite);
      } else {
        new_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
      }
    } else {
      person->determine_clinical_or_not(new_parasite);
    }

  }

  person->schedule_mature_gametocyte_event(new_parasite);
}
