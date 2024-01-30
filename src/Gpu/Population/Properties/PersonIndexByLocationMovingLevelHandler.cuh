/* 
 * File:   PersonIndexByLocationMovingLevelHandler.h
 * Author: Merlin
 *
 * Created on August 1, 2013, 9:27 PM
 */

#ifndef PERSONINDEXBYLOCATIONMOVINGLEVELHANDLER_CUH
#define    PERSONINDEXBYLOCATIONMOVINGLEVELHANDLER_CUH

#include "IndexHandler.cuh"
#include "Core/PropertyMacro.h"

namespace GPU {
    class PersonIndexByLocationMovingLevelHandler;
}

class GPU::PersonIndexByLocationMovingLevelHandler : public GPU::IndexHandler {
 public:
  PersonIndexByLocationMovingLevelHandler();

  //    PersonIndexByLocationMovingLevelHandler(const PersonIndexByLocationMovingLevelHandler& orig);
  virtual ~PersonIndexByLocationMovingLevelHandler();

 private:

};

#endif    /* PERSONINDEXBYLOCATIONMOVINGLEVELHANDLER_H */

