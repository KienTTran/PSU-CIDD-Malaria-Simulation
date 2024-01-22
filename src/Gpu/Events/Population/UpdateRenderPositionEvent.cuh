//
// Created by kient on 8/4/2023.
//

#ifndef PERSONUPDATERENDERPOSITIONEVENT_CUH
#define PERSONUPDATERENDERPOSITIONEVENT_CUH


#include "Gpu/Events/Event.cuh"

namespace GPU {
    class UpdateRenderPositionEvent;
}

class GPU::UpdateRenderPositionEvent: public GPU::Event {

DISALLOW_COPY_AND_ASSIGN(UpdateRenderPositionEvent)

public:
    UpdateRenderPositionEvent();

    virtual ~UpdateRenderPositionEvent();

    static void schedule_event(GPU::Scheduler *scheduler, const int &time);

    std::string name() override;

private:
    void execute() override;


};


#endif //PERSONUPDATERENDERPOSITIONEVENT_CUH
