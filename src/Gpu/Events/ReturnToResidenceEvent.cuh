/* 
 * File:   ReturnToResidenceEvent.h
 * Author: Merlin
 *
 * Created on August 2, 2013, 11:20 AM
 */

#ifndef RETURNTORESIDENCEEVENT_CUH
#define    RETURNTORESIDENCEEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include <string>

namespace GPU {
    class Person;
    class Scheduler;
    class ReturnToResidenceEvent;
}

class GPU::ReturnToResidenceEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ReturnToResidenceEvent)

 DISALLOW_MOVE(ReturnToResidenceEvent)


 public:
  ReturnToResidenceEvent();

  //    ReturnToResidenceEvent(const ReturnToResidenceEvent& orig);
  virtual ~ReturnToResidenceEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &time);

  std::string name() override {
    return "ReturnToResidenceEvent";
  }

 private:
  void execute() override;
};

#endif    /* RETURNTORESIDENCEEVENT_CUH */
