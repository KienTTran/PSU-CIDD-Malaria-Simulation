/* 
 * File:   Dispatcher.cpp
 * Author: nguyentran
 * 
 * Created on May 3, 2013, 3:46 PM
 */

#include "Dispatcher.cuh"
#include "Gpu/Population/Properties/IndexHandler.cuh"
#include "Gpu/Events/Event.cuh"
#include "Helpers/ObjectHelpers.h"

GPU::Dispatcher::Dispatcher() : events_(nullptr) {}

void GPU::Dispatcher::init() {
  events_ = new GPUEventPtrVector();
}

GPU::Dispatcher::~Dispatcher() {
  GPU::Dispatcher::clear_events();
  ObjectHelpers::delete_pointer<GPUEventPtrVector>(events_);
}

void GPU::Dispatcher::add(GPU::Event *event) {
  events_->push_back(event);
  event->IndexHandler::set_index(events_->size() - 1);
}

void GPU::Dispatcher::remove(GPU::Event *event) {
  events_->back()->IndexHandler::set_index(event->IndexHandler::index());

  //    std::cout << "1"<<event->IndexHandler::index()<< std::endl;
  events_->at(event->IndexHandler::index()) = events_->back();
  //    std::cout << "2"<< std::endl;

  events_->pop_back();
  event->IndexHandler::set_index(-1);
}

void GPU::Dispatcher::clear_events() {
  if (events_==nullptr) return;
  if (events_->empty()) return;
  //    std::cout << "Clear event"<< std::endl;

  for (auto &event : *events_) {
    event->dispatcher = nullptr;
    event->executable = false;
  }

  events_->clear();
}

void GPU::Dispatcher::update() {
  events_->shrink_to_fit();
}
