/*
 * File:   Scheduler.cpp
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:27 PM
 */

#include "Scheduler.cuh"

#include <vector>

#include "Core/Config/Config.h"
#include "Dispatcher.cuh"
#include "Helpers/ObjectHelpers.h"
#include "Helpers/TimeHelpers.h"
#include "Model.h"
#include "easylogging++.h"
#include "Gpu/Events/Event.cuh"
#include "Gpu/Events/Population/UpdateByLocationEvent.cuh"
#include "Gpu/Events/Population/UpdateRenderPositionEvent.cuh"
#include "Gpu/Events/Population/UpdateRenderOGLEvent.cuh"

using namespace date;

GPU::Scheduler::Scheduler(Model* model)
    : current_time_(-1), total_available_time_(-1), model_(model), is_force_stop_(false) {}

GPU::Scheduler::~Scheduler() {
  clear_all_events();
}

void GPU::Scheduler::extend_total_time(int new_total_time) {
  if (total_available_time_ < new_total_time) {
    for (auto i = total_available_time_; i <= new_total_time; i++) {
      individual_events_list_.push_back(GPUEventPtrVector());
      population_events_list_.push_back(GPUEventPtrVector());
    }
  }
  total_available_time_ = new_total_time;
}

void GPU::Scheduler::clear_all_events() {
  clear_all_events(individual_events_list_);
  clear_all_events(population_events_list_);
}

void GPU::Scheduler::initialize(const date::year_month_day& starting_date, const int& total_time) {
  // 720 here is to prevent schedule birthday event at the end of simulation
  set_total_available_time(total_time + 720);
  set_current_time(0);

  calendar_date = sys_days(starting_date);
}

void GPU::Scheduler::clear_all_events(GPUEventPtrVector2& events_list) const {
  for (auto& timestep_events : events_list) {
    for (auto* event : timestep_events) {
      if (event->dispatcher != nullptr) {
        event->dispatcher->remove(event);
      }
      ObjectHelpers::delete_pointer<Event>(event);
    }
    timestep_events.clear();
  }
  events_list.clear();
}

int GPU::Scheduler::total_available_time() const {
  return total_available_time_;
}

void GPU::Scheduler::set_total_available_time(const int& value) {
  if (total_available_time_ > 0) {
    clear_all_events();
  }
  total_available_time_ = value;
  individual_events_list_.assign(total_available_time_, GPUEventPtrVector());
  population_events_list_.assign(total_available_time_, GPUEventPtrVector());
}

void GPU::Scheduler::schedule_individual_event(GPU::Event* event) {
  schedule_event(individual_events_list_[event->time], event);
}

void GPU::Scheduler::schedule_population_event(GPU::Event* event) {
  if (event->time < population_events_list_.size()) {
    schedule_event(population_events_list_[event->time], event);
  }
}

void GPU::Scheduler::schedule_event(GPUEventPtrVector& time_events, Event* event) {
  // Schedule event in the future
  // Event time cannot exceed total time or less than current time
  if (event->time > Model::CONFIG->total_time() || event->time < current_time_) {
    LOG_IF(event->time < current_time_, FATAL)
        << "Error when schedule event " << event->name() << " at " << event->time << ". Current_time: " << current_time_
        << " - total time: " << total_available_time_;
    VLOG(2) << "Cannot schedule event " << event->name() << " at " << event->time << ". Current_time: " << current_time_
            << " - total time: " << total_available_time_;
    ObjectHelpers::delete_pointer<Event>(event);
  } else {
    time_events.push_back(event);
    event->scheduler = this;
    event->executable = true;
  }
}

void GPU::Scheduler::cancel(GPU::Event* event) {
  event->executable = false;
}

void GPU::Scheduler::execute_events_list(GPUEventPtrVector& events_list) const {
  for (auto& event : events_list) {
    // std::cout << event->name() << std::endl;
    event->perform_execute();
    ObjectHelpers::delete_pointer<Event>(event);
  }
  ObjectHelpers::clear_vector_memory<Event>(events_list);
}

void GPU::Scheduler::run() {
  LOG(INFO) << "Simulation is running";
  current_time_ = 0;

  for (current_time_ = 0; !can_stop(); current_time_++) {
    LOG_IF(current_time_ % Model::CONFIG->debug_config().log_interval == 0, INFO) << "Day: " << current_time_;
    begin_time_step();

    daily_update();

    end_time_step();
    calendar_date += days { 1 };
  }
  Model::MODEL->model_finished = true;
}

void GPU::Scheduler::begin_time_step() const {
  if (model_ != nullptr) {
    model_->begin_time_step();
  }
}

void GPU::Scheduler::daily_update() {
  if (model_ != nullptr) {
    model_->daily_update();

    if (is_today_first_day_of_month()) {
      GPU::UpdateByLocationEvent::schedule_event(Model::GPU_SCHEDULER, current_time_);
      model_->monthly_update();
    }

    if (is_today_first_day_of_year()) {
      // std::cout << date::year_month_day{calendar_date} << std::endl;
      model_->yearly_update();
    }

    if (Model::CONFIG->render_config().display_gui) {
        GPU::UpdateRenderPositionEvent::schedule_event(Model::GPU_SCHEDULER, current_time_);
        GPU::UpdateRenderOGLEvent::schedule_event(Model::GPU_SCHEDULER, current_time_);
    }

    // population related events
    execute_events_list(population_events_list_[current_time_]);

    // individual related events
    execute_events_list(individual_events_list_[current_time_]);
  }
}

void GPU::Scheduler::end_time_step() const {
  if (model_ != nullptr) {
    model_->end_time_step();
  }
}

bool GPU::Scheduler::can_stop() const {
  return current_time_ > Model::CONFIG->total_time() || is_force_stop_;
}

int GPU::Scheduler::current_day_in_year() const {
  return TimeHelpers::day_of_year(calendar_date);
}

bool GPU::Scheduler::is_today_last_day_of_year() const {
  year_month_day ymd { calendar_date };
  return ymd.month() == month { 12 } && ymd.day() == day { 31 };
}

bool GPU::Scheduler::is_today_first_day_of_month() const {
  // return true;
  year_month_day ymd { calendar_date };
  return ymd.day() == day { 1 };
}

bool GPU::Scheduler::is_today_first_day_of_year() const {
  year_month_day ymd { calendar_date };
  return ymd.month() == month { 1 } && ymd.day() == day { 1 };
}

bool GPU::Scheduler::is_today_last_day_of_month() const {
  const auto next_date = calendar_date + days { 1 };
  year_month_day ymd { next_date };
  return ymd.day() == day { 1 };
}
