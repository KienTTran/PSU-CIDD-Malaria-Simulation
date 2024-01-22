/* 
 * File:   PersonIndexAllHandler.h
 * Author: nguyentran
 *
 * Created on April 17, 2013, 10:21 AM
 */

#ifndef PERSONINDEXALLHANDLER_CUH
#define    PERSONINDEXALLHANDLER_CUH

#include "Core/PropertyMacro.h"
#include "IndexHandler.cuh"

namespace GPU {
    class PersonIndexAllHandler;
}
class GPU::PersonIndexAllHandler : public GPU::IndexHandler {
 public:
  PersonIndexAllHandler();

  virtual ~PersonIndexAllHandler();

 private:

};

#endif    /* PERSONINDEXALLHANDLER_H */

