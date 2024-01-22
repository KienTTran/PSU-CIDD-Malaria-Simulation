/* 
 * File:   ImportationEvent.h
 * Author: Merlin
 *
 * Created on February 21, 2014, 2:42 PM
 */

#ifndef IMPORTATIONPERIODICALLYEVENT_CUH
#define    IMPORTATIONPERIODICALLYEVENT_CUH

#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"

namespace GPU {
    class ImportationPeriodicallyEvent;
}

class GPU::ImportationPeriodicallyEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ImportationPeriodicallyEvent)


 VIRTUAL_PROPERTY_REF(int, location)

 VIRTUAL_PROPERTY_REF(int, duration)

 VIRTUAL_PROPERTY_REF(int, genotype_id)

 VIRTUAL_PROPERTY_REF(int, number_of_cases)

 public:
  ImportationPeriodicallyEvent(const int &location = -1, const int &duration = -1, int genotype_id = -1,
                               const int &number_of_cases = -1, const int &start_day = -1);

  //    ImportationEvent(const ImportationEvent& orig);
  virtual ~ImportationPeriodicallyEvent();

  static void schedule_event(GPU::Scheduler *scheduler, const int &location, const int &duration, unsigned int genotype_id,
                             const int &number_of_cases, const int &start_day);

  std::string name() override {
    return "ImportationPeriodicallyEvent";
  }

 private:
  void execute() override;

};

#endif    /* IMPORTATIONPERIODICALLYEVENT_CUH */
