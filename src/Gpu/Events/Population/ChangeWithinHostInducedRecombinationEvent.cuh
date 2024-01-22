//
// Created by kient on 5/3/2022.
//

#ifndef POMS_CHANGEWITHINHOSTINDUCEDRECOMBINATIONEVENT_CUH
#define POMS_CHANGEWITHINHOSTINDUCEDRECOMBINATIONEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include <string>
#include <vector>
#include <easylogging++.h>
#include "Model.h"
#include "Core/Config/Config.h"

namespace GPU{
    class ChangeWithinHostInducedRecombinationEvent;
}

class GPU::ChangeWithinHostInducedRecombinationEvent : public GPU::Event {
  DISALLOW_COPY_AND_ASSIGN(ChangeWithinHostInducedRecombinationEvent)
  DISALLOW_MOVE(ChangeWithinHostInducedRecombinationEvent)
public:
  bool value {true};
public:
  explicit ChangeWithinHostInducedRecombinationEvent(const bool &value = 0.0, const int &at_time = -1);

  ~ChangeWithinHostInducedRecombinationEvent() override = default;

  std::string name() override {
    return "ChangeWithinHostInducedRecombinationEvent";
  }

private:
  void execute() override;
};

#endif  // POMS_CHANGEWITHINHOSTINDUCEDRECOMBINATIONEVENT_CUH
