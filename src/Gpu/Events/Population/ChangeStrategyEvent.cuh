#ifndef CHANGESTRATEGYEVENT_CUH
#define CHANGESTRATEGYEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU {
    class ChangeStrategyEvent;
}

class GPU::ChangeStrategyEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ChangeStrategyEvent)

 DISALLOW_MOVE(ChangeStrategyEvent)

 public:
  int strategy_id{-1};

  ChangeStrategyEvent(const int &at_time, const int &strategy_id);

  virtual ~ChangeStrategyEvent() = default;

  std::string name() override {
    return "ChangeStrategyEvent";
  }

 private:
  void execute() override;
};

#endif // CHANGESTRATEGYEVENT_CUH
