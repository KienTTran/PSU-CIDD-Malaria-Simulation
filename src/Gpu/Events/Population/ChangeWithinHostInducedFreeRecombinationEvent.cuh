//
// Created by kient on 5/3/2022.
//

#ifndef POMS_CHANGEWITHINHOSTINDUCEDFREERECOMBINATIONEVENT_CUH
#define POMS_CHANGEWITHINHOSTINDUCEDFREERECOMBINATIONEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include <string>
#include <vector>
#include <easylogging++.h>
#include "Model.h"
#include "Core/Config/Config.h"

namespace GPU {
    class ChangeWithinHostInducedFreeRecombinationEvent;
}

class GPU::ChangeWithinHostInducedFreeRecombinationEvent : public GPU::Event {
  DISALLOW_COPY_AND_ASSIGN(ChangeWithinHostInducedFreeRecombinationEvent)
  DISALLOW_MOVE(ChangeWithinHostInducedFreeRecombinationEvent)
public:
  bool value {true};
public:
  explicit ChangeWithinHostInducedFreeRecombinationEvent(const bool &value = false, const int &at_time = -1);

  ~ChangeWithinHostInducedFreeRecombinationEvent() override = default;

  std::string name() override {
    return "ChangeWithinHostInducedRecombinationEvent";
  }

private:
  void execute() override;
};

#endif  // POMS_CHANGEWITHINHOSTINDUCEDFREERECOMBINATIONEVENT_H
