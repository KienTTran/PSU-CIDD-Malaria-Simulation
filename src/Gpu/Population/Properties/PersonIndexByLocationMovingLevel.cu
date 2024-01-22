/* 
* File:   PersonIndexByLocationMovingLevel.cpp
* Author: Merlin
* 
* Created on August 1, 2013, 9:26 PM
*/

#include "PersonIndexByLocationMovingLevel.cuh"
#include "Core/Config/Config.h"
#include "Model.h"
#include <cassert>

GPU::PersonIndexByLocationMovingLevel::PersonIndexByLocationMovingLevel(const int &no_location, const int &no_level) {
  Initialize(no_location, no_level);
}

GPU::PersonIndexByLocationMovingLevel::~PersonIndexByLocationMovingLevel() {
  vPerson_.clear();
}

void GPU::PersonIndexByLocationMovingLevel::Initialize(const int &no_location, const int &no_level) {
  vPerson_.clear();

  GPUPersonPtrVector ppv;
  GPUPersonPtrVector2 ppv2;
  ppv2.assign(no_level, ppv);

  vPerson_.assign(no_location, ppv2);

  Model::CONFIG->h_popsize_by_moving_level = IntVector(no_location*no_level,0);
}

void GPU::PersonIndexByLocationMovingLevel::add(GPU::Person *p) {
  assert(vPerson_.size() > p->location() && p->location() >= 0);
  assert(vPerson_[p->location()].size() > p->moving_level());
  add(p, p->location(), p->moving_level());
}

void GPU::PersonIndexByLocationMovingLevel::remove(GPU::Person *p) {
  remove_without_set_index(p);
  p->GPU::PersonIndexByLocationMovingLevelHandler::set_index(-1);
}

void GPU::PersonIndexByLocationMovingLevel::notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue,
                                                     const void *newValue) {
  switch (property) {
    case GPU::Person::LOCATION:change_property(p, *(int *) newValue, p->moving_level());
      break;
    case GPU::Person::MOVING_LEVEL:change_property(p, p->location(), *(int *) newValue);
      break;
    default:break;
  }
}

std::size_t GPU::PersonIndexByLocationMovingLevel::size() const {
  return 0;
}

void GPU::PersonIndexByLocationMovingLevel::add(GPU::Person *p, const int &location, const int &moving_level) {
  vPerson_[location][moving_level].push_back(p);
  p->GPU::PersonIndexByLocationMovingLevelHandler::set_index(vPerson_[location][moving_level].size() - 1);
  Model::CONFIG->h_popsize_by_moving_level[location * Model::CONFIG->circulation_info().number_of_moving_levels + moving_level] = vPerson_[location][moving_level].size();
}

void GPU::PersonIndexByLocationMovingLevel::remove_without_set_index(GPU::Person *p) {
  vPerson_[p->location()][p->moving_level()].back()->GPU::PersonIndexByLocationMovingLevelHandler::set_index(
      p->GPU::PersonIndexByLocationMovingLevelHandler::index());
  vPerson_[p->location()][p->moving_level()][p->GPU::PersonIndexByLocationMovingLevelHandler::index()] =
      vPerson_[p->location()][p->moving_level()].back();
  vPerson_[p->location()][p->moving_level()].pop_back();
  Model::CONFIG->h_popsize_by_moving_level[p->location() * Model::CONFIG->circulation_info().number_of_moving_levels + p->moving_level()] = vPerson_[p->location()][p->moving_level()].size();
}

void GPU::PersonIndexByLocationMovingLevel::change_property(GPU::Person *p, const int &location, const int &moving_level) {
  //remove from old position
  remove_without_set_index(p); //to save 1 set and improve performance since the index of p will changed when add

  //add to new position
  add(p, location, moving_level);
}

void GPU::PersonIndexByLocationMovingLevel::update() {
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    for (int ml = 0; ml < Model::CONFIG->circulation_info().number_of_moving_levels; ml++) {
      std::vector<Person *>(vPerson_[location][ml]).swap(vPerson_[location][ml]);
    }
  }
}
