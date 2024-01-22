/* 
 * File:   UpdateByHavingDrugEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 12:16 PM
 */

#ifndef UPDATEWHENDRUGISPRESENTEVENT_CUH
#define    UPDATEWHENDRUGISPRESENTEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include <string>

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class UpdateWhenDrugIsPresentEvent;
}



class GPU::UpdateWhenDrugIsPresentEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(UpdateWhenDrugIsPresentEvent)


 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

 public:
  UpdateWhenDrugIsPresentEvent();

  //    UpdateByHavingDrugEvent(const UpdateByHavingDrugEvent& orig);
  virtual ~UpdateWhenDrugIsPresentEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                             const int &time);

  std::string name() override {
    return "UpdateByHavingDrugEvent";
  }

 private:
  void execute() override;
};

#endif    /* UPDATEWHENDRUGISPRESENTEVENT_CUH */
