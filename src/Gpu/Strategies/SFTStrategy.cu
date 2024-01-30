/* 
 * File:   SFTStrategy.cpp
 * Author: nguyentran
 * 
 * Created on June 3, 2013, 8:00 PM
 */

#include <cassert>
#include "SFTStrategy.cuh"
#include "Gpu/Therapies/Therapy.cuh"
#include "IStrategy.cuh"
#include <sstream>

GPU::SFTStrategy::SFTStrategy() : GPU::IStrategy("SFTStrategy", SFT) {}

GPU::SFTStrategy::~SFTStrategy() = default;

std::vector<GPU::Therapy *> &GPU::SFTStrategy::get_therapy_list() {
  return therapy_list_;
}

void GPU::SFTStrategy::add_therapy(GPU::Therapy *therapy) {
  therapy_list_.push_back(therapy);
}

GPU::Therapy *GPU::SFTStrategy::get_therapy(GPU::Person *person) {
  return therapy_list_[0];
}

std::string GPU::SFTStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << "-" << therapy_list_[0]->id();
  return sstm.str();
}

void GPU::SFTStrategy::update_end_of_time_step() {
  //do nothing here
}

void GPU::SFTStrategy::adjust_started_time_point(const int &current_time) {
  //do nothing
}

void GPU::SFTStrategy::monthly_update() {
  //do nothing
}
