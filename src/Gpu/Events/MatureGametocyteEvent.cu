/* 
 * File:   MatureGametocyteEvent.cu
 * Author: Merlin
 * 
 * Created on July 31, 2013, 11:38 PM
 */

#include "MatureGametocyteEvent.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Population/SingleHostClonalParasitePopulations.h"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Person.cuh"

GPU::MatureGametocyteEvent::MatureGametocyteEvent() : blood_parasite_(nullptr) {}

GPU::MatureGametocyteEvent::~MatureGametocyteEvent() = default;

void GPU::MatureGametocyteEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *blood_parasite,
                                           const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new GPU::MatureGametocyteEvent();
    e->dispatcher = p;
    e->set_blood_parasite(blood_parasite);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::MatureGametocyteEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  LOG_IF(person->index() >= 1000 && person->index() <= 1085,INFO)
    << fmt::format("GPU::MatureGametocyteEvent::execute() {}",person->index());
  if (person->all_clonal_parasite_populations()->contain(blood_parasite_)) {
    blood_parasite_->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
  }

}
