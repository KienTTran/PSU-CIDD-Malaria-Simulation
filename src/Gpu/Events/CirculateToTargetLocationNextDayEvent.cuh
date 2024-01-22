/* 
 * File:   CirculateToTargetLocationNextDayEvent.h
 * Author: Merlin
 *
 * Created on August 2, 2013, 10:57 AM
 */

#ifndef CIRCULATETOTARGETLOCATIONNEXTDAYEVENT_CUH
#define    CIRCULATETOTARGETLOCATIONNEXTDAYEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"

class Scheduler;

namespace GPU{
    class Person;
    class CirculateToTargetLocationNextDayEvent;
}

class GPU::CirculateToTargetLocationNextDayEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(CirculateToTargetLocationNextDayEvent)

 DISALLOW_MOVE(CirculateToTargetLocationNextDayEvent)


 PROPERTY_REF(int, target_location)

 public:
  CirculateToTargetLocationNextDayEvent();

  //    CirculateToTargetLocationNextDayEvent(const CirculateToTargetLocationNextDayEvent& orig);
  virtual ~CirculateToTargetLocationNextDayEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &target_location, const int &time);

  std::string name() override;

 private:
  void execute() override;

 private:

};

#endif    /* CIRCULATETOTARGETLOCATIONNEXTDAYEVENT_CUH */
