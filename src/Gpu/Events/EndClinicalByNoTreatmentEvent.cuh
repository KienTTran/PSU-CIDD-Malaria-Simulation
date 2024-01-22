/* 
 * File:   EndClinicalByNoTreatmentEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 12:28 PM
 */

#ifndef ENDCLINICALBYNOTREATMENTEVENT_CUH
#define    ENDCLINICALBYNOTREATMENTEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class EndClinicalByNoTreatmentEvent;
    class Scheduler;
}

class GPU::EndClinicalByNoTreatmentEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(EndClinicalByNoTreatmentEvent)


 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

 public:
  EndClinicalByNoTreatmentEvent();

  //    EndClinicalByNoTreatmentEvent(const EndClinicalByNoTreatmentEvent& orig);
  virtual ~EndClinicalByNoTreatmentEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                             const int &time);

  std::string name() override {
    return "EndClinicalByNoTreatmentEvent";
  }

 private:
  void execute() override;
};

#endif    /* ENDCLINICALBYNOTREATMENTEVENT_CUH */
