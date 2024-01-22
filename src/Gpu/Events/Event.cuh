/* 
 * File:   Event.h
 * Author: nguyentran
 *
 * Created on May 3, 2013, 3:13 PM
 */

#ifndef EVENT_CUH
#define    EVENT_CUH

#include "Gpu/Population/Properties/IndexHandler.cuh"
#include "Core/PropertyMacro.h"
#include <string>


namespace GPU{
    class Event;
    class Dispatcher;
    class Scheduler;
}

class GPU::Event : public GPU::IndexHandler {
 public:
  GPU::Scheduler *scheduler{nullptr};
  GPU::Dispatcher *dispatcher{nullptr};
  bool executable{false};
  int time{-1};

  Event();

  //    Event(const Event& orig);
  virtual ~Event();

  void perform_execute();

  virtual std::string name() = 0;

 private:
  virtual void execute() = 0;

};

#endif    /* EVENT_H */
