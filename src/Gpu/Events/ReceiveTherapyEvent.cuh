/* 
 * File:   ReceiveTherapyEvent.h
 * Author: Merlin
 *
 * Created on November 4, 2014, 3:00 PM
 */

#ifndef RECEIVETHERAPYEVENT_CUH
#define    RECEIVETHERAPYEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Population/ClonalParasitePopulation.h"



class Therapy;

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class ReceiveTherapyEvent;
}

class GPU::ReceiveTherapyEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ReceiveTherapyEvent)

 DISALLOW_MOVE(ReceiveTherapyEvent)

 POINTER_PROPERTY(Therapy, received_therapy)

 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

public:
  bool is_part_of_MAC_therapy { false };

public:
  ReceiveTherapyEvent();

  //    ReceiveTherapyEvent(const ReceiveTherapyEvent& orig);
  virtual ~ReceiveTherapyEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, Therapy *therapy, const int &time,
                             GPU::ClonalParasitePopulation *clinical_caused_parasite, bool is_part_of_MAC_therapy = false);

  std::string name() override {
    return "ReceiveTherapyEvent";
  }

 private:
  void execute() override;
};

#endif    /* RECEIVETHERAPYEVENT_CUH */
