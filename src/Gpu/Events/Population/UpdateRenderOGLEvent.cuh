//
// Created by kient on 8/4/2023.
//

#ifndef MASS_PERSONUPDATERENDEROGLEVENT_CUH
#define MASS_PERSONUPDATERENDEROGLEVENT_CUH


#include "Events/Event.h"

class UpdateRenderOGLEvent: public Event {

DISALLOW_COPY_AND_ASSIGN(UpdateRenderOGLEvent)

public:
    UpdateRenderOGLEvent();

    virtual ~UpdateRenderOGLEvent();

    static void schedule_event(Scheduler *scheduler, const int &time);

    std::string name() override;

private:
    void execute() override;

};


#endif //MASS_PERSONUPDATERENDEROGLEVENT_CUH
