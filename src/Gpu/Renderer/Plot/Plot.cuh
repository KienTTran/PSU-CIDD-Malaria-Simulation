//
// Created by ktt on 2/20/24.
//

#ifndef POMS_PLOT_CUH
#define POMS_PLOT_CUH
#include <iostream>
#include "Core/TypeDef.h"
#include "Core/PropertyMacro.h"
#include <imgui.h>
#include "implot.h"
#include "implot_internal.h"

namespace GPU{
    class Plot;
}

class GPU::Plot {
public:
    PROPERTY_REF(TVector<GPU::Plot*>,plots);
public:
    Plot();
    ~Plot();
    void init();
    void add_plot(GPU::Plot* plot);
    virtual void start_plot();
    virtual void update_plot();
};


#endif //POMS_PLOT_CUH
