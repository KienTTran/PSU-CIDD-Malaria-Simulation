/*
 * BurkinaFaso.hxx
 * 
 * This events class contains events that needs to be run for Burkina Faso at 
 * the end of the calendar year.
 */
#ifndef BURKINAFASOEVENTS_HXX
#define BURKINAFASOEVENTS_HXX

#include <math.h>

#include "easylogging++.h"

#include "Core/Config/Config.h"
#include "Core/Scheduler.h"
#include "Events/Event.h"
#include "Helpers/TimeHelpers.h"
#include "Model.h"

class AnnualCoverageUpdateEvent : public Event {
  private:
    float rate_ = 0.0;

    // Execute the annual events for Burkina Faso.
    void execute() override {

      // Grab a reference to the location_db to work with
      auto& location_db = Model::CONFIG->location_db();

      // Iterate though and adjust the rates
      auto count = Model::CONFIG->number_of_locations();
      for (auto ndx = 0; ndx < count; ndx++) {
        location_db[ndx].p_treatment_less_than_5 = adjust(location_db[ndx].p_treatment_less_than_5, rate_);
        location_db[ndx].p_treatment_more_than_5 = adjust(location_db[ndx].p_treatment_more_than_5, rate_);
      }

      // Schedule for one year from now
      auto time = scheduler->current_time() + TimeHelpers::number_of_days_to_next_year(scheduler->calendar_date);
      AnnualCoverageUpdateEvent* event = new AnnualCoverageUpdateEvent(rate_, time);
      scheduler->schedule_population_event(event);

      // Log on demand
      VLOG(1) << "Updated coverage by " << rate_ << ", scheduled next for: " << time;
    }

    // Update the coverage by the given rate, round up at the third decimal place but clamped at 1.0 (100%)
    double adjust(double coverage, double rate) {
      coverage += (1 - coverage) * rate;
      coverage = int(coverage * pow(10, 3)) / pow(10, 3);
      return (coverage < 1.0) ? coverage : 1.0;
    }

  public:
    inline static const std::string EventName = "annual_coverage_update_event";

    AnnualCoverageUpdateEvent(float rate, int start) : rate_(rate) { time = start; }
    ~AnnualCoverageUpdateEvent() = default;

    // Return the name of this event
    std::string name() override { return EventName; }

};

#endif