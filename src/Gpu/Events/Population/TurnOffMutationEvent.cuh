#ifndef TURNOFFMUTATIONEVENT_CUH
#define TURNOFFMUTATIONEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU {
    class TurnOffMutationEvent;
}

class GPU::TurnOffMutationEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(TurnOffMutationEvent)

 DISALLOW_MOVE(TurnOffMutationEvent)

 public:
  explicit TurnOffMutationEvent(const int &at_time);

  virtual ~TurnOffMutationEvent() = default;

  std::string name() override {
    return "TurnOffMutationEvent";
  }

 private:
  void execute() override;
};

#endif // TURNOFFMUTATIONEVENT_CUH
