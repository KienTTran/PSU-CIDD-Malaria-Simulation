/* 
 * File:   SwitchImmuneComponentEvent.cpp
 * Author: Merlin
 * 
 * Created on July 3, 2013, 3:06 PM
 */

#include "SwitchImmuneComponentEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Population/NonInfantImmuneComponent.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/Person.cuh"
#include <cassert>

GPU::SwitchImmuneComponentEvent::SwitchImmuneComponentEvent() = default;

GPU::SwitchImmuneComponentEvent::~SwitchImmuneComponentEvent() = default;

void GPU::SwitchImmuneComponentEvent::execute() {

  assert(dispatcher!=nullptr);
  auto *p = dynamic_cast<GPU::Person *>(dispatcher);
  p->immune_system()->set_immune_component(new GPU::NonInfantImmuneComponent());

}

void GPU::SwitchImmuneComponentEvent::schedule_for_switch_immune_component_event(Scheduler *scheduler, GPU::Person *p,
                                                                            const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new SwitchImmuneComponentEvent();
    e->dispatcher = p;
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}
