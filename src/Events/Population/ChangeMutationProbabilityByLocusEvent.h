//
// Created by kient on 8/22/2023.
//

#ifndef POMS_CHANGEMUTATIONPROBABILITYBYLOCUSEVENT_H
#define POMS_CHANGEMUTATIONPROBABILITYBYLOCUSEVENT_H

#include "Core/PropertyMacro.h"
#include "Events/Event.h"
#include "Core/Scheduler.h"
#include <string>
#include <vector>
#include <easylogging++.h>
#include "Model.h"
#include "Core/Config/Config.h"

class ChangeMutationProbabilityByLocusEvent : public Event {
    DISALLOW_COPY_AND_ASSIGN(ChangeMutationProbabilityByLocusEvent)
    DISALLOW_MOVE(ChangeMutationProbabilityByLocusEvent)
public:
    double value {0.001};
public:
    explicit ChangeMutationProbabilityByLocusEvent(const double &value = 0.0, const int &at_time = -1);

    ~ChangeMutationProbabilityByLocusEvent() override = default;

    std::string name() override {
        return "ChangeMutationProbabilityByLocusEvent";
    }

private:
    void execute() override;
};


#endif //POMS_CHANGEMUTATIONPROBABILITYBYLOCUSEVENT_H
