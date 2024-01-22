/* 
 * File:   BirthdayEvent.h
 * Author: nguyentran
 *
 * Created on May 9, 2013, 2:42 PM
 */

#ifndef BIRTHDAYEVENT_CUH
#define    BIRTHDAYEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Core/ObjectPool.h"
#include <string>

namespace GPU{
    class Person;
    class BirthdayEvent;
}

class GPU::BirthdayEvent : public GPU::Event {

 DISALLOW_COPY_AND_ASSIGN(BirthdayEvent)

 public:
  BirthdayEvent();

  //    BirthdayEvent(const BirthdayEvent& orig);
  virtual ~BirthdayEvent();

  static void schedule_event(GPU::Scheduler *scheduler, GPU::Person *p, const int &time);

  std::string name() override;

 private:
  void execute() override;
};

#endif    /* BIRTHDAYEVENT_CUH */
