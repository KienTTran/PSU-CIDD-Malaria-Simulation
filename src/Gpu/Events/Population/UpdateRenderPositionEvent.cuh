//
// Created by kient on 8/4/2023.
//

#ifndef MASS_PERSONUPDATERENDERPOSITIONEVENT_CUH
#define MASS_PERSONUPDATERENDERPOSITIONEVENT_CUH


#include "Events/Event.h"

class UpdateRenderPositionEvent: public Event {

DISALLOW_COPY_AND_ASSIGN(UpdateRenderPositionEvent)

public:
    UpdateRenderPositionEvent();

    virtual ~UpdateRenderPositionEvent();

    static void schedule_event(Scheduler *scheduler, const int &time);

    std::string name() override;

private:
    void execute() override;


};


#endif //MASS_PERSONUPDATERENDERPOSITIONEVENT_CUH
