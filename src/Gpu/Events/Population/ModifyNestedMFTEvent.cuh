#ifndef MODIFYNESTEDMFTEVENT_CUH
#define MODIFYNESTEDMFTEVENT_CUH

#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU {
    class ModifyNestedMFTEvent;
}

class GPU::ModifyNestedMFTEvent : public GPU::Event {
 DISALLOW_COPY_AND_ASSIGN(ModifyNestedMFTEvent)

 DISALLOW_MOVE(ModifyNestedMFTEvent)

 public:
  int strategy_id{-1};

  ModifyNestedMFTEvent(const int &at_time, const int &strategy_id);

  virtual ~ModifyNestedMFTEvent() = default;

  std::string name() override {
    return "ChangeStrategyEvent";
  }

 private:
  void execute() override;
};

#endif // MODIFYNESTEDMFTEVENT_CUH
