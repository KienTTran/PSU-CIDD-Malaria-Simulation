#include "ModifyNestedMFTEvent.cuh"
#include "Model.h"
#include "Gpu/Strategies/IStrategy.cuh"
#include "Gpu/Strategies/NestedMFTMultiLocationStrategy.cuh"
#include "Core/Config/Config.h"
#include "Gpu/Strategies/NestedMFTStrategy.cuh"

GPU::ModifyNestedMFTEvent::ModifyNestedMFTEvent(const int &at_time, const int &strategy_id) : strategy_id(strategy_id) {
  time = at_time;
}

void GPU::ModifyNestedMFTEvent::execute() {
  GPU::IStrategy *new_strategy = nullptr;
  if (Model::GPU_TREATMENT_STRATEGY->get_type()==GPU::IStrategy::NestedMFTMultiLocation) {
    new_strategy = Model::CONFIG->gpu_strategy_db()[strategy_id];
    dynamic_cast<GPU::NestedMFTMultiLocationStrategy *>(Model::GPU_TREATMENT_STRATEGY)->strategy_list[0] = new_strategy;
//    new_strategy->adjust_started_time_point(Model::GPU_SCHEDULER->current_time());
  }

  if (Model::GPU_TREATMENT_STRATEGY->get_type()==GPU::IStrategy::NestedMFT) {
    new_strategy = Model::CONFIG->gpu_strategy_db()[strategy_id];
    dynamic_cast<GPU::NestedMFTStrategy *>(Model::GPU_TREATMENT_STRATEGY)->strategy_list[0] = new_strategy;
    new_strategy->adjust_started_time_point(Model::GPU_SCHEDULER->current_time());
  }

  if (new_strategy==nullptr) {
    LOG(FATAL) << "Modify Nested MFT Event error with null ptr.";
  }

  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : switch first strategy in nested strategies to "
            << new_strategy->name;
}
