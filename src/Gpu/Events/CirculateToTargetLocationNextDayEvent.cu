/* 
 * File:   CirculateToTargetLocationNextDayEvent.cpp
 * Author: Merlin
 * 
 * Created on August 2, 2013, 10:57 AM
 */

#include "CirculateToTargetLocationNextDayEvent.cuh"
#include "Model.h"
#include "Core/Random.h"
#include "Core/Config/Config.h"
#include "ReturnToResidenceEvent.cuh"
#include "Gpu/Population/Person.cuh"


GPU::CirculateToTargetLocationNextDayEvent::CirculateToTargetLocationNextDayEvent() : target_location_(0) {}

GPU::CirculateToTargetLocationNextDayEvent::~CirculateToTargetLocationNextDayEvent() = default;

void GPU::CirculateToTargetLocationNextDayEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &target_location,
                                                           const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new CirculateToTargetLocationNextDayEvent();
    e->dispatcher = p;
    e->set_target_location(target_location);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

std::string GPU::CirculateToTargetLocationNextDayEvent::name() {
  return "CirculateToTargetLocationNextDayEvent";
}

void GPU::CirculateToTargetLocationNextDayEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  if(person->index() >= 1040 && person->index() <= 1045){
      printf("GPU::CirculateToTargetLocationNextDayEvent::execute() %d\n",person->index());
  }
  person->set_location(target_location_);

  if (target_location_!=person->residence_location()) {
    //if person already have return trip then no need to reschedule it
    //elase
    //Schedule for a return trip in next several days base on gamma distribution
    if (!person->has_return_to_residence_event()) {
      int length_of_trip = 0;
      while (length_of_trip < 1) {
        length_of_trip = static_cast<int>(std::round(
            Model::RANDOM->random_gamma(Model::CONFIG->circulation_info().length_of_stay_theta,
                                        Model::CONFIG->circulation_info().length_of_stay_k)));
      }

      //            std::cout << length_of_trip << std::endl;
      GPU::ReturnToResidenceEvent::schedule_event(Model::GPU_SCHEDULER, person,
                                             Model::GPU_SCHEDULER->current_time() + length_of_trip);

    }
  } else {
    //return by chance so we cancel all return event
    //cancel return trip and do nothing
    person->cancel_all_return_to_residence_events();
  }

}
