#pragma once

#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU {
    class IntroduceLumefantrineMutantEvent;
}

class GPU::IntroduceLumefantrineMutantEvent : public GPU::Event {
  DISALLOW_COPY_AND_ASSIGN(IntroduceLumefantrineMutantEvent)

    DISALLOW_MOVE(IntroduceLumefantrineMutantEvent)


    VIRTUAL_PROPERTY_REF(int, location)

    VIRTUAL_PROPERTY_REF(double, fraction)

public:
  explicit IntroduceLumefantrineMutantEvent(const int& location = -1, const int& execute_at = -1,
    const double& fraction = 0);

  //    ImportationEvent(const ImportationEvent& orig);
  ~IntroduceLumefantrineMutantEvent() override;

  std::string name() override {
    return "IntroduceLumefantrineMutantEvent";
  }

private:
  void execute() override;

};
