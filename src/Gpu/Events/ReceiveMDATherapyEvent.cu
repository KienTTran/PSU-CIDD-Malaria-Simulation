/* 
 * File:   ReceiveMDADrugEvent.cpp
 * Author: Merlin
 * 
 * Created on August 25, 2013, 10:20 PM
 */
#include "ReceiveMDATherapyEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Therapies/Therapy.h"
#include "Model.h"
#include "Gpu/Population/Person.cuh"

GPU::ReceiveMDATherapyEvent::ReceiveMDATherapyEvent() : received_therapy_(nullptr) {};

GPU::ReceiveMDATherapyEvent::~ReceiveMDATherapyEvent() = default;

void GPU::ReceiveMDATherapyEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, Therapy *therapy, const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new ReceiveMDATherapyEvent();
    e->dispatcher = p;
    e->set_received_therapy(therapy);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::ReceiveMDATherapyEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  //    if (person->is_in_external_population()) {
  //        return;
  //    }

  person->receive_therapy(received_therapy_, nullptr);

  //if this person has progress to clinical event then cancel it
  person->cancel_all_other_progress_to_clinical_events_except(nullptr);
  person->change_all_parasite_update_function(Model::MODEL->gpu_progress_to_clinical_update_function(),
                                              Model::MODEL->gpu_immunity_clearance_update_function());

  person->schedule_update_by_drug_event(nullptr);

}
