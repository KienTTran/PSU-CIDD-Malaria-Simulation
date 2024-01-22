/* 
 * File:   Scheduler.h
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:27 PM
 */

#ifndef SCHEDULER_CUH
#define  SCHEDULER_CUH

#include <chrono>
#include "date/date.h"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

namespace GPU{
    class Event;
    class Scheduler;
}

class Model;

class GPU::Scheduler {
 DISALLOW_COPY_AND_ASSIGN(Scheduler)

 DISALLOW_MOVE(Scheduler)

 PROPERTY_REF(int, current_time)

 PROPERTY_HEADER(int, total_available_time)

 POINTER_PROPERTY(Model, model)

 PROPERTY_REF(bool, is_force_stop)

 public:
  date::sys_days calendar_date;

  GPUEventPtrVector2 individual_events_list_;
  GPUEventPtrVector2 population_events_list_;

  explicit Scheduler(Model *model = nullptr);

  virtual ~Scheduler();

  void extend_total_time(int new_total_time);

  void clear_all_events();

  void clear_all_events(GPUEventPtrVector2 &events_list) const;

  virtual void schedule_individual_event(GPU::Event *event);

  virtual void schedule_population_event(GPU::Event *event);

  virtual void schedule_event(GPUEventPtrVector &time_events, Event *event);

  virtual void cancel(GPU::Event *event);

  void execute_events_list(GPUEventPtrVector &events_list) const;

  void initialize(const date::year_month_day &starting_date, const int &total_time);

  void run();

  void begin_time_step() const;

  void end_time_step() const;

  bool can_stop() const;

  int current_day_in_year() const;

  bool is_today_last_day_of_month() const;

  bool is_today_first_day_of_month() const;

  bool is_today_first_day_of_year() const;

  bool is_today_last_day_of_year() const;

  void daily_update();
};

#endif  /* SCHEDULER_H */
