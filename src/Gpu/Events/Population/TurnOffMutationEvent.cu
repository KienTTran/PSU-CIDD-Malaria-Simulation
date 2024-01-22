#include "TurnOffMutationEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "easylogging++.h"
#include "Model.h"
#include "Core/Config/Config.h"

GPU::TurnOffMutationEvent::TurnOffMutationEvent(const int &at_time) {
  time = at_time;
}

void GPU::TurnOffMutationEvent::execute() {
  Model::CONFIG->mutation_probability_by_locus() = 0.0;
  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : turn mutation off";
}
