#include "TurnOnMutationEvent.cuh"
#include "easylogging++.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

GPU::TurnOnMutationEvent::TurnOnMutationEvent(const int &at_time, const double &mutation_probability) :
  mutation_probability{mutation_probability} {
  time = at_time;
}

void GPU::TurnOnMutationEvent::execute() {
  Model::CONFIG->mutation_probability_by_locus() = mutation_probability;
  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : turn mutation on";
}
