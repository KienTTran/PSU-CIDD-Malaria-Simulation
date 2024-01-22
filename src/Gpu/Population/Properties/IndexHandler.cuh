/* 
 * File:   PersonIndexHandler.h
 * Author: nguyentran
 *
 * Created on April 17, 2013, 10:29 AM
 */

#ifndef PERSONINDEXHANDLER_CUH
#define    PERSONINDEXHANDLER_CUH

#include "Core/PropertyMacro.h"

namespace GPU {
    class IndexHandler;
}

class GPU::IndexHandler {
 DISALLOW_COPY_AND_ASSIGN(IndexHandler)

 PROPERTY_REF(std::size_t, index)

 public:
  IndexHandler();

  virtual ~IndexHandler();

};

#endif    /* PERSONINDEXHANDLER_H */

