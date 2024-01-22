/* 
 * File:   EndClinicalEvent.cpp
 * Author: Merlin
 * 
 * Created on July 31, 2013, 12:27 PM
 */

#include "EndClinicalEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/Person.cuh"

GPU::EndClinicalEvent::EndClinicalEvent() : clinical_caused_parasite_(nullptr) {}

GPU::EndClinicalEvent::~EndClinicalEvent() = default;

void
GPU::EndClinicalEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                                 const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new EndClinicalEvent();
    e->dispatcher = p;
    e->set_clinical_caused_parasite(clinical_caused_parasite);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::EndClinicalEvent::execute() {
  auto person = dynamic_cast<GPU::Person *>(dispatcher);

  if (person->all_clonal_parasite_populations()->size()==0) {
    person->change_state_when_no_parasite_in_blood();

  } else {
    //still have parasite in blood
    person->immune_system()->set_increase(true);
    person->set_host_state(GPU::Person::ASYMPTOMATIC);

    person->determine_relapse_or_not(clinical_caused_parasite_);

  }
}
