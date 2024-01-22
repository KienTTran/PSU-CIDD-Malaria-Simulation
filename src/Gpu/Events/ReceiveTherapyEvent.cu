/* 
 * File:   ReceiveTherapyEvent.cpp
 * Author: Merlin
 * 
 * Created on November 4, 2014, 3:00 PM
 */

#include "ReceiveTherapyEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Therapies/Therapy.h"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/Person.cuh"

GPU::ReceiveTherapyEvent::ReceiveTherapyEvent() : received_therapy_(nullptr), clinical_caused_parasite_(nullptr) {}

GPU::ReceiveTherapyEvent::~ReceiveTherapyEvent() = default;

void GPU::ReceiveTherapyEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, Therapy *therapy, const int &time,
                                         GPU::ClonalParasitePopulation *clinical_caused_parasite,
                                         bool is_part_of_MAC_therapy) {
  if (scheduler != nullptr) {
    auto *e = new ReceiveTherapyEvent();
    e->dispatcher = p;
    e->set_received_therapy(therapy);
    e->time = time;
    e->set_clinical_caused_parasite(clinical_caused_parasite);
    e->is_part_of_MAC_therapy = is_part_of_MAC_therapy;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::ReceiveTherapyEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  //    if (person->is_in_external_population()) {
  //        return;
  //    }

  person->receive_therapy(received_therapy_, clinical_caused_parasite_, is_part_of_MAC_therapy);

  person->schedule_update_by_drug_event(clinical_caused_parasite_);
}
