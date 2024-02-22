/* 
 * File:   ReturnToResidenceEvent.cu
 * Author: Merlin
 * 
 * Created on August 2, 2013, 11:20 AM
 */

#include <cassert>

#include "ReturnToResidenceEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Population/Person.cuh"


GPU::ReturnToResidenceEvent::ReturnToResidenceEvent() = default;

GPU::ReturnToResidenceEvent::~ReturnToResidenceEvent() = default;

void GPU::ReturnToResidenceEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new GPU::ReturnToResidenceEvent();
    e->dispatcher = p;
    e->time = time;
    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::ReturnToResidenceEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  person->set_location(person->residence_location());

}
