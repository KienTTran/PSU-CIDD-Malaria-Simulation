/* 
 * File:   ProgressToClinicalEvent.h
 * Author: Merlin
 *
 * Created on July 30, 2013, 2:36 PM
 */

#ifndef PROGRESSTOCLINICALEVENT_CUH
#define    PROGRESSTOCLINICALEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include <string>

namespace GPU{
    class ClonalParasitePopulation;
    class Person;
    class Scheduler;
    class ProgressToClinicalEvent;
}

class GPU::ProgressToClinicalEvent : public GPU::Event {

 DISALLOW_COPY_AND_ASSIGN(ProgressToClinicalEvent)

 DISALLOW_MOVE(ProgressToClinicalEvent)

 POINTER_PROPERTY(GPU::ClonalParasitePopulation, clinical_caused_parasite)

 public:
  ProgressToClinicalEvent();

  virtual ~ProgressToClinicalEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                             const int &time);

  static void receive_no_treatment_routine(GPU::Person *p);

  std::string name() override {
    return "ProgressToClinicalEvent";
  }

 private:
  void execute() override;
};

#endif    /* PROGRESSTOCLINICALEVENT_CUH */
