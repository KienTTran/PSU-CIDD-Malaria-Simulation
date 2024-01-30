/* 
 * File:   MFTStrategy.cpp
 * Author: nguyentran
 * 
 * Created on June 4, 2013, 11:09 AM
 */

#include "MFTStrategy.cuh"
#include "Core/Random.h"
#include "Model.h"
#include <sstream>
#include "IStrategy.cuh"
#include "Gpu/Therapies/Therapy.cuh"

GPU::MFTStrategy::MFTStrategy() : GPU::IStrategy("MFTStrategy", MFT) {}

GPU::MFTStrategy::~MFTStrategy() = default;

void GPU::MFTStrategy::add_therapy(GPU::Therapy *therapy) {
  therapy_list.push_back(therapy);
}

GPU::Therapy *GPU::MFTStrategy::get_therapy(GPU::Person *person) {

  const auto p = Model::RANDOM->random_flat(0.0, 1.0);

  double sum = 0;
  for (auto i = 0; i < distribution.size(); i++) {
    sum += distribution[i];
    if (p <= sum) {
      return therapy_list[i];
    }
  }

  return therapy_list[therapy_list.size() - 1];
}

std::string GPU::MFTStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << "-";
  std::string sep;
  for (auto *therapy : therapy_list) {
    sstm << sep << therapy->id();
    sep = ",";
  }
  sep = "";
  sstm << "-";
  for (auto dist : distribution) {
    sstm << sep << dist;
    sep = ",";
  }
  return sstm.str();
}

void GPU::MFTStrategy::adjust_started_time_point(const int &current_time) {}

void GPU::MFTStrategy::monthly_update() {
  //do nothing here

}

void GPU::MFTStrategy::update_end_of_time_step() {
  //do nothing here
}
