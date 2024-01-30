/* 
 * File:   CyclingStrategy.cpp
 * Author: nguyentran
 * 
 * Created on June 4, 2013, 11:10 AM
 */

#include "CyclingStrategy.cuh"
#include "Model.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Core/Config/Config.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include <sstream>
#include "IStrategy.cuh"
#include "Gpu/Therapies/Therapy.cuh"

GPU::CyclingStrategy::CyclingStrategy() : GPU::IStrategy("CyclingStrategy", Cycling) {}

GPU::CyclingStrategy::~CyclingStrategy() = default;

void GPU::CyclingStrategy::add_therapy(GPU::Therapy *therapy) {
  therapy_list.push_back(therapy);
}

void GPU::CyclingStrategy::switch_therapy() {
  //    std::cout << "Switch from: " << index_ << "\t - to: " << index_ + 1;
  index++;
  index %= therapy_list.size();
  Model::GPU_DATA_COLLECTOR->update_UTL_vector();

  // TODO: cycling_time should be match with calendar day
  next_switching_day = Model::GPU_SCHEDULER->current_time() + cycling_time;
  LOG(INFO) << date::year_month_day{Model::GPU_SCHEDULER->calendar_date}
            << ": Cycling Strategy switch therapy to: " << therapy_list[index]->id();
}

GPU::Therapy *GPU::CyclingStrategy::get_therapy(GPU::Person *person) {

  //int index = ((Global::scheduler->currentTime - Global::startTreatmentDay) / circleTime) % therapyList.size();
  //    std::cout << therapy_list()[index_]->id() << std::endl;
  return therapy_list[index];
}

std::string GPU::CyclingStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << "-";
  std::string sep;
  for (auto *therapy : therapy_list) {
    sstm << sep << therapy->id();
    sep = ",";
  }
  return sstm.str();
}

void GPU::CyclingStrategy::update_end_of_time_step() {
  if (Model::GPU_SCHEDULER->current_time()==next_switching_day) {
    switch_therapy();
    //            std::cout << to_string() << std::endl;
  }
}

void GPU::CyclingStrategy::adjust_started_time_point(const int &current_time) {
  next_switching_day = Model::GPU_SCHEDULER->current_time() + cycling_time;
  index = 0;
}

void GPU::CyclingStrategy::monthly_update() {}
