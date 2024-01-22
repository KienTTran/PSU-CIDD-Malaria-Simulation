/* 
 * File:   MoveParasiteToBloodEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 11:14 PM
 */

#ifndef MOVEPARASITETOBLOODEVENT_CUH
#define    MOVEPARASITETOBLOODEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include <string>

namespace GPU{
    class ClonalParasitePopulation;
    class Genotype;
    class Person;
    class Scheduler;
    class MoveParasiteToBloodEvent;
}

class GPU::MoveParasiteToBloodEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(MoveParasiteToBloodEvent)

 DISALLOW_MOVE(MoveParasiteToBloodEvent)


 POINTER_PROPERTY(GPU::Genotype, infection_genotype)

 public:
  MoveParasiteToBloodEvent();

  //    MoveParasiteToBloodEvent(const MoveParasiteToBloodEvent& orig);
  virtual ~MoveParasiteToBloodEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::Genotype *infection_type, const int &time);

  std::string name() override {
    return "MoveParasiteToBloodEvent";
  }

 private:
  void execute() override;
};

#endif    /* MOVEPARASITETOBLOODEVENT_CUH */
