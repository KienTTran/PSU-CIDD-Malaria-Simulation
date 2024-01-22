/* 
 * File:   PersonIndexByLocationStateAgeClassHandler.h
 * Author: nguyentran
 *
 * Created on May 2, 2013, 10:57 AM
 */

#ifndef PERSONINDEXBYLOCATIONSTATEAGECLASSHANDLER_CUH
#define    PERSONINDEXBYLOCATIONSTATEAGECLASSHANDLER_CUH

#include "Core/PropertyMacro.h"
#include "IndexHandler.cuh"

namespace GPU {
    class PersonIndexByLocationStateAgeClassHandler;
}

class GPU::PersonIndexByLocationStateAgeClassHandler : public GPU::IndexHandler {

 public:
  PersonIndexByLocationStateAgeClassHandler();

//    PersonIndexByLocationStateAgeClassHandler(const PersonIndexByLocationStateAgeClassHandler& orig);
  virtual ~PersonIndexByLocationStateAgeClassHandler();

 private:

};

#endif    /* PERSONINDEXBYLOCATIONSTATEAGECLASSHANDLER_H */

