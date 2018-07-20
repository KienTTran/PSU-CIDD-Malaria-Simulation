/* 
 * File:   ReturnToResidenceEvent.cpp
 * Author: Merlin
 * 
 * Created on August 2, 2013, 11:20 AM
 */

#include <cassert>

#include "ReturnToResidenceEvent.h"
#include "Person.h"
#include "Scheduler.h"

OBJECTPOOL_IMPL(ReturnToResidenceEvent)

ReturnToResidenceEvent::ReturnToResidenceEvent() {
}

ReturnToResidenceEvent::~ReturnToResidenceEvent() {
}

void ReturnToResidenceEvent::schedule_event(Scheduler* scheduler, Person* p, const int& time) {
    if (scheduler != nullptr) {
        ReturnToResidenceEvent* e = new ReturnToResidenceEvent();
        e->dispatcher = p;

        e->executable = true;
        e->time = time;

        
        p->add(e);
        scheduler->schedule_individual_event(e);
    }
}

void ReturnToResidenceEvent::execute() {
    Person* person = (Person*) dispatcher;
    person->set_location(person->residence_location());

}
