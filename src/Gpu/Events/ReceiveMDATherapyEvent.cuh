/* 
 * File:   ReceiveMDADrugEvent.h
 * Author: Merlin
 *
 * Created on August 25, 2013, 10:20 PM
 */

#ifndef RECEIVEMDADRUGEVENT_CUH
#define    RECEIVEMDADRUGEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include <string>
#include "Core/PropertyMacro.h"


namespace GPU {
    class Person;
    class Scheduler;
    class ReceiveMDATherapyEvent;
    class Therapy;
}


class GPU::ReceiveMDATherapyEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ReceiveMDATherapyEvent)

 DISALLOW_MOVE(ReceiveMDATherapyEvent)

 POINTER_PROPERTY(GPU::Therapy, received_therapy)

 public:
  ReceiveMDATherapyEvent();

  //    ReceiveMDADrugEvent(const ReceiveMDADrugEvent& orig);
  virtual ~ReceiveMDATherapyEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::Therapy *therapy, const int &time);

  std::string name() override {
    return "ReceiveMDADrugEvent";
  }

 private:
  void execute() override;

};

#endif    /* RECEIVEMDADRUGEVENT_CUH */
