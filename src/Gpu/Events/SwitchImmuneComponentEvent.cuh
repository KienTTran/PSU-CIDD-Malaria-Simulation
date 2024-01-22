/* 
 * File:   SwitchImmuneComponentEvent.h
 * Author: Merlin
 *
 * Created on July 3, 2013, 3:06 PM
 */

#ifndef SWITCHIMMUNECOMPONENTEVENT_CUH
#define    SWITCHIMMUNECOMPONENTEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"


namespace GPU {
    class Person;
    class Scheduler;
    class SwitchImmuneComponentEvent;
}

class GPU::SwitchImmuneComponentEvent : public GPU::Event {

 public:
  SwitchImmuneComponentEvent();

  SwitchImmuneComponentEvent(const SwitchImmuneComponentEvent &orig);

  virtual ~SwitchImmuneComponentEvent();

  static void schedule_for_switch_immune_component_event(Scheduler *scheduler, GPU::Person *p, const int &time);

  virtual std::string name() {
    return "SwitchImmuneComponentEvent";
  }

 private:
  virtual void execute();
};

#endif    /* SWITCHIMMUNECOMPONENTEVENT_CUH */

