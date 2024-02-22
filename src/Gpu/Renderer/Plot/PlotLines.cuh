//
// Created by ktt on 2/20/24.
//

#ifndef POMS_PLOTLINES_CUH
#define POMS_PLOTLINES_CUH
#include "Plot.cuh"

namespace GPU{
    class PlotLines;
}

// utility structure for realtime plot
struct ScrollingBuffer {
    int MaxSize;
    int Offset;
    ImVector<ImVec2> Data;
    ScrollingBuffer(int max_size = 2000) {
      MaxSize = max_size;
      Offset  = 0;
      Data.reserve(MaxSize);
    }
    void AddPoint(float x, float y) {
      if (Data.size() < MaxSize)
        Data.push_back(ImVec2(x,y));
      else {
        Data[Offset] = ImVec2(x,y);
        Offset =  (Offset + 1) % MaxSize;
      }
    }
    void Erase() {
      if (Data.size() > 0) {
        Data.shrink(0);
        Offset  = 0;
      }
    }
};

class GPU::PlotLines : public GPU::Plot{
public:
    PlotLines();
    ~PlotLines();
    void start_plot() override;
    void update_plot() override;
private:
    ImPlotRect limits;
    ScrollingBuffer sdata1;
    float t = 0;
    double *xs1, *ys1, *ys2, *ys3;
};


#endif //POMS_PLOTLINES_CUH
