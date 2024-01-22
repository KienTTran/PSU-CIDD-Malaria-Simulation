/* 
 * File:   PersonIndexAll.cpp
 * Author: nguyentran
 * 
 * Created on April 17, 2013, 10:15 AM
 */

#include <vector>
#include "PersonIndexAll.cuh"

GPU::PersonIndexAll::PersonIndexAll() = default;

GPU::PersonIndexAll::~PersonIndexAll() {
  vPerson_.clear();

}

void GPU::PersonIndexAll::add(GPU::Person *p) {
  vPerson_.push_back(p);
  p->GPU::PersonIndexAllHandler::set_index(vPerson_.size() - 1);
}

void GPU::PersonIndexAll::remove(GPU::Person *p) {
  //move the last element to current position and remove the last holder
  vPerson_.back()->GPU::PersonIndexAllHandler::set_index(p->GPU::PersonIndexAllHandler::index());
  vPerson_[p->GPU::PersonIndexAllHandler::index()] = vPerson_.back();
  vPerson_.pop_back();
  p->GPU::PersonIndexAllHandler::set_index(-1);
  //    delete p;
  //    p = nullptr;
}

std::size_t GPU::PersonIndexAll::size() const {
  return vPerson_.size();
}

void GPU::PersonIndexAll::notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue,
                                   const void *newValue) {}

void GPU::PersonIndexAll::update() {
  vPerson_.shrink_to_fit();
}
