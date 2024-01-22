/* 
 * File:   MatureGametocyteEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 11:38 PM
 */

#ifndef MATUREGAMETOCYTEEVENT_CUH
#define    MATUREGAMETOCYTEEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class MatureGametocyteEvent;
}

class GPU::MatureGametocyteEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(MatureGametocyteEvent);

 POINTER_PROPERTY(GPU::ClonalParasitePopulation, blood_parasite)

 public:
  MatureGametocyteEvent();

  //    MatureGametocyteEvent(const MatureGametocyteEvent& orig);
  virtual ~MatureGametocyteEvent();

  static void
  schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *blood_parasite, const int &time);

  std::string name() override {
    return "MatureGametocyteEvent";
  }

 private:
  void execute() override;
};

#endif    /* MATUREGAMETOCYTEEVENT_CUH */
