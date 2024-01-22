/* 
 * File:   Dispatcher.h
 * Author: nguyentran
 *
 * Created on May 3, 2013, 3:46 PM
 */

#ifndef DISPATCHER_CUH
#define    DISPATCHER_CUH

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

namespace GPU{
    class Event;
    class Dispatcher;
}

class GPU::Dispatcher {
 DISALLOW_COPY_AND_ASSIGN(Dispatcher)

 POINTER_PROPERTY(GPUEventPtrVector, events)

 public:
  Dispatcher();

  //    Dispatcher(const Dispatcher& orig);
  virtual ~Dispatcher();

  virtual void init();

  virtual void add(GPU::Event *event);

  virtual void remove(GPU::Event *event);

  virtual void clear_events();

  virtual void update();

};

#endif    /* DISPATCHER_H */

