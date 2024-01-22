/* 
 * File:   EndClinicalDueToDrugResistanceEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 11:24 AM
 */

#ifndef ENDCLINICALDUETODRUGRESISTANCEEVENT_CUH
#define    ENDCLINICALDUETODRUGRESISTANCEEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class EndClinicalDueToDrugResistanceEvent;
}

class GPU::EndClinicalDueToDrugResistanceEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(EndClinicalDueToDrugResistanceEvent)


 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

 public:
  EndClinicalDueToDrugResistanceEvent();

  //    EndClinicalDueToDrugResistanceEvent(const EndClinicalDueToDrugResistanceEvent& orig);
  virtual ~EndClinicalDueToDrugResistanceEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                             const int &time);

  std::string name() override {
    return "EndClinicalDueToDrugResistanceEvent";
  }

 private:
  void execute() override;

};

#endif    /* ENDCLINICALDUETODRUGRESISTANCEEVENT_CUH */
