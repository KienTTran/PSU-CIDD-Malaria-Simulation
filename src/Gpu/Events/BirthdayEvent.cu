/* 
 * File:   BirthdayEvent.cpp
 * Author: nguyentran
 * 
 * Created on May 9, 2013, 2:42 PM
 */

#include <cassert>

#include "BirthdayEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Helpers/TimeHelpers.h"
#include "easylogging++.h"
#include "Gpu/Population/Person.cuh"

GPU::BirthdayEvent::BirthdayEvent() = default;

GPU::BirthdayEvent::~BirthdayEvent() = default;

void GPU::BirthdayEvent::execute() {
  assert(dispatcher!=nullptr);
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  person->increase_age_by_1_year();

  const auto days_to_next_year = TimeHelpers::number_of_days_to_next_year(scheduler->calendar_date);

  schedule_event(scheduler, person, scheduler->current_time() + days_to_next_year);
}

void GPU::BirthdayEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &time) {
  if (scheduler!=nullptr) {
    auto *birthday_event = new BirthdayEvent();
    birthday_event->dispatcher = p;
    birthday_event->time = time;

    p->add(birthday_event);
    scheduler->schedule_individual_event(birthday_event);
    //        std::cout << scheduler->current_time() << " - hello" << std::endl;
  }
}

std::string GPU::BirthdayEvent::name() {
  return "Birthday Event";
}
