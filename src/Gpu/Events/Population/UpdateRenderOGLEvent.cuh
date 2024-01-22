//
// Created by kient on 8/4/2023.
//

#ifndef PERSONUPDATERENDEROGLEVENT_CUH
#define PERSONUPDATERENDEROGLEVENT_CUH

#include "Gpu/Events/Event.cuh"

namespace GPU {
    class UpdateRenderOGLEvent;
}

class GPU::UpdateRenderOGLEvent: public GPU::Event {

DISALLOW_COPY_AND_ASSIGN(UpdateRenderOGLEvent)

public:
    UpdateRenderOGLEvent();

    virtual ~UpdateRenderOGLEvent();

    static void schedule_event(GPU::Scheduler *scheduler, const int &time);

    std::string name() override;

private:
    void execute() override;

};


#endif //PERSONUPDATERENDEROGLEVENT_CUH
