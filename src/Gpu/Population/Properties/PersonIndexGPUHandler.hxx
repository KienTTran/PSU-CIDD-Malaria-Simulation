//
// Created by kient on 12/9/2023.
//

#ifndef PERSONINDEXGPUHANDLER_HXX
#define PERSONINDEXGPUHANDLER_HXX

#include "Core/PropertyMacro.h"
#include "IndexHandler.cuh"

namespace GPU {
    class PersonIndexGPUHandler;
}

class GPU::PersonIndexGPUHandler : public GPU::IndexHandler {
DISALLOW_COPY_AND_ASSIGN(PersonIndexGPUHandler)

public:
    PersonIndexGPUHandler(){}

    virtual ~PersonIndexGPUHandler(){}

private:

};


#endif //PERSONINDEXGPUHANDLER_HXX
