//
// Created by kient on 12/24/2023.
//

#ifndef PERSONUPDATEBYLOCATIONEVENT_CUH
#define PERSONUPDATEBYLOCATIONEVENT_CUH

#include "Gpu/Events/Event.cuh"

namespace GPU {
    class UpdateByLocationEvent;
}

class GPU::UpdateByLocationEvent : public GPU::Event {

DISALLOW_COPY_AND_ASSIGN(UpdateByLocationEvent)

public:
    UpdateByLocationEvent();

virtual ~UpdateByLocationEvent();

static void schedule_event(GPU::Scheduler *scheduler, const int &time);

std::string name() override;

private:
    void execute() override;


};


#endif //PERSONUPDATEBYLOCATIONEVENT_CUH
