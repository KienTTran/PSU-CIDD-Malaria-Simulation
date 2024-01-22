#ifndef TURNONMUTATIONEVENT_CUH
#define TURNONMUTATIONEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU {
    class TurnOnMutationEvent;
}

class GPU::TurnOnMutationEvent : public GPU::Event {
DISALLOW_COPY_AND_ASSIGN(TurnOnMutationEvent)

DISALLOW_MOVE(TurnOnMutationEvent)

  double mutation_probability = 0.0;
  int drug_id = -1;

public:
  explicit TurnOnMutationEvent(const int &at_time, const double &mutation_probability);

  ~TurnOnMutationEvent() override = default;

  std::string name() override {
    return "TurnOnMutationEvent";
  }

private:
  void execute() override;
};

#endif // TURNONMUTATIONEVENT_CUH
