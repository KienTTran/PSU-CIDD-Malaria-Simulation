#include "ChangeStrategyEvent.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Strategies/IStrategy.cuh"

GPU::ChangeStrategyEvent::ChangeStrategyEvent(const int &at_time, const int &strategy_id) : strategy_id(strategy_id) {
  time = at_time;
}

void GPU::ChangeStrategyEvent::execute() {
  Model::MODEL->set_treatment_strategy(strategy_id);
  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << ": switch to " << Model::GPU_TREATMENT_STRATEGY->name;
}
