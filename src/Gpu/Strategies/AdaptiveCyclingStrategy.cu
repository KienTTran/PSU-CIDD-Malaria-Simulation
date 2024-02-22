/* 
 * File:   AdaptiveCyclingStrategy.cu
 * Author: nguyentran
 * 
 * Created on June 4, 2013, 11:10 AM
 */

#include "AdaptiveCyclingStrategy.cuh"
#include "Model.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Core/Config/Config.h"
#include <sstream>
#include "IStrategy.cuh"
#include "Gpu/Therapies/Therapy.cuh"

GPU::AdaptiveCyclingStrategy::AdaptiveCyclingStrategy() : GPU::IStrategy("AdaptiveCyclingStrategy", AdaptiveCycling) {}

GPU::AdaptiveCyclingStrategy::~AdaptiveCyclingStrategy() = default;

void GPU::AdaptiveCyclingStrategy::add_therapy(GPU::Therapy *therapy) {
  therapy_list.push_back(therapy);
}

void GPU::AdaptiveCyclingStrategy::switch_therapy() {
  //    std::cout << "Switch from: " << index_ << "\t - to: " << index_ + 1;
  index++;
  index %= therapy_list.size();

  Model::GPU_DATA_COLLECTOR->update_UTL_vector();
  LOG(INFO) << date::year_month_day{Model::GPU_SCHEDULER->calendar_date}
            << ": Adaptive Cycling Strategy switch Therapy to: " << therapy_list[index]->id();
}

GPU::Therapy *GPU::AdaptiveCyclingStrategy::get_therapy(GPU::Person *person) {
  return therapy_list[index];
}

std::string GPU::AdaptiveCyclingStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << "-";
  std::string sep;
  for (auto *therapy : therapy_list) {
    sstm << sep << therapy->id();
    sep = ",";
  }
  return sstm.str();
}

void GPU::AdaptiveCyclingStrategy::update_end_of_time_step() {

  if (Model::GPU_SCHEDULER->current_time()==latest_switch_time) {
    switch_therapy();
    //            std::cout << to_string() << std::endl;
  } else {
    if (Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()[get_therapy(nullptr)->id()] > trigger_value) {
      // TODO:: turn_off_days and delay_until_actual_trigger should be match with calendar day
      if (Model::GPU_SCHEDULER->current_time() > latest_switch_time + turn_off_days) {
        latest_switch_time = Model::GPU_SCHEDULER->current_time() + delay_until_actual_trigger;
        LOG(INFO) << date::year_month_day{Model::GPU_SCHEDULER->calendar_date}
                  << ": Adaptive Cycling will switch therapy next year";
        //                    std::cout << "TF: " << Model::GPU_DATA_COLLECTOR->current_TF_by_therapy()[get_therapy()->id()] << std::endl;
      }
    }
  }

}

void GPU::AdaptiveCyclingStrategy::adjust_started_time_point(const int &current_time) {
  latest_switch_time = -1;
  index = 0;
}

void GPU::AdaptiveCyclingStrategy::monthly_update() {}
