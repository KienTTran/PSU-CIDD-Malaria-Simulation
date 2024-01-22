//
// Created by Nguyen Tran on 2019-05-01.
//

#ifndef PCMS_INTRODUCEPLAS2COPYPARASITEEVENT_CUH
#define PCMS_INTRODUCEPLAS2COPYPARASITEEVENT_CUH


#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Gpu/Events/Event.cuh"
#include <string>

namespace GPU{
    class IntroducePlas2CopyParasiteEvent;
}

class GPU::IntroducePlas2CopyParasiteEvent : public GPU::Event {
DISALLOW_COPY_AND_ASSIGN(IntroducePlas2CopyParasiteEvent)

DISALLOW_MOVE(IntroducePlas2CopyParasiteEvent)


VIRTUAL_PROPERTY_REF(int, location)

VIRTUAL_PROPERTY_REF(double, fraction)

public:
  explicit IntroducePlas2CopyParasiteEvent(const int &location = -1, const int &execute_at = -1,
                                           const double &fraction = 0);

  //    ImportationEvent(const ImportationEvent& orig);
  ~IntroducePlas2CopyParasiteEvent() override;

  std::string name() override {
    return "Plas2ImportationEvent";
  }

private:
  void execute() override;

};

#endif //PCMS_INTRODUCEPLAS2COPYPARASITEEVENT_CUH
