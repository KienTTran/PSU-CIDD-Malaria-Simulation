/* 
 * File:   EndClinicalEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 12:27 PM
 */

#ifndef ENDCLINICALEVENT_CUH
#define    ENDCLINICALEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class EndClinicalEvent;
}

class GPU::EndClinicalEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(EndClinicalEvent);

 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

 public:
  EndClinicalEvent();

  //    EndClinicalEvent(const EndClinicalEvent& orig);
  virtual ~EndClinicalEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                             const int &time);

  std::string name() override {
    return "EndClinicalEvent";
  }

 private:
  void execute() override;
};

#endif    /* ENDCLINICALEVENT_CUH */
