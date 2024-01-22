//
// Created by nguyentd on 11/3/2020.
//

#ifndef POMS_SRC_EVENTS_POPULATION_INTRODUCETRIPLEMUTANTTODPMEVENT_CUH
#define POMS_SRC_EVENTS_POPULATION_INTRODUCETRIPLEMUTANTTODPMEVENT_CUH


#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU{
    class IntroduceTrippleMutantToDPMEvent;
}

class GPU::IntroduceTrippleMutantToDPMEvent : public GPU::Event {
DISALLOW_COPY_AND_ASSIGN(IntroduceTrippleMutantToDPMEvent)

DISALLOW_MOVE(IntroduceTrippleMutantToDPMEvent)


VIRTUAL_PROPERTY_REF(int, location)

VIRTUAL_PROPERTY_REF(double, fraction)

public:
  explicit IntroduceTrippleMutantToDPMEvent(
      const int& location = -1, const int& execute_at = -1,
      const double& fraction = 0
  );

  //    ImportationEvent(const ImportationEvent& orig);
  ~IntroduceTrippleMutantToDPMEvent() override;

  std::string name() override {
    return "IntroduceTrippleMutantToDPMEvent";
  }

private:
  void execute() override;

};


#endif //POMS_SRC_EVENTS_POPULATION_INTRODUCETRIPLEMUTANTTODPMEVENT_CUH
