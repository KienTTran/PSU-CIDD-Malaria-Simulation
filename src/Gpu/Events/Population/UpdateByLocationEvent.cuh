//
// Created by kient on 12/24/2023.
//

#ifndef MASS_PERSONUPDATEBYLOCATIONEVENT_CUH
#define MASS_PERSONUPDATEBYLOCATIONEVENT_CUH

#include "Events/Event.h"

class UpdateByLocationEvent : public Event {

DISALLOW_COPY_AND_ASSIGN(UpdateByLocationEvent)

public:
    UpdateByLocationEvent();

virtual ~UpdateByLocationEvent();

static void schedule_event(Scheduler *scheduler, const int &time);

std::string name() override;

private:
    void execute() override;


};


#endif //MASS_PERSONUPDATEBYLOCATIONEVENT_CUH
