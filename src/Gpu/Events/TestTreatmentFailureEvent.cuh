/* 
 * File:   TestTreatmentFailureEvent.h
 * Author: Merlin
 *
 * Created on July 31, 2013, 11:36 AM
 */

#ifndef TESTTREATMENTFAILUREEVENT_CUH
#define    TESTTREATMENTFAILUREEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class TestTreatmentFailureEvent;
}



class GPU::TestTreatmentFailureEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(TestTreatmentFailureEvent)


 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)
  //    PROPERTY_REF(bool, isResistance);
 PROPERTY_REF(int, therapyId)

 public:
  TestTreatmentFailureEvent();

  //    TestTreatmentFailureEvent(const TestTreatmentFailureEvent& orig);
  virtual ~TestTreatmentFailureEvent();

  static void
  schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite, const int &time,
                 const int &t_id = 0);

  std::string name() override {
    return "TestTreatmentFailureEvent";
  }

 private:
  void execute() override;

};

#endif    /* TESTTREATMENTFAILUREEVENT_CUH */
