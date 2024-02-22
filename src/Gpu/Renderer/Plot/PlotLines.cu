//
// Created by ktt on 2/20/24.
//
#include "PlotLines.cuh"
#include <imgui_stdlib.h>
#include "Model.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Core/Config/Config.h"

GPU::PlotLines::~PlotLines() {

}

GPU::PlotLines::PlotLines() {

}

void GPU::PlotLines::start_plot() {
  xs1 = new double[Model::CONFIG->total_time()];
  ys1 = new double[Model::CONFIG->total_time()];
  ys2 = new double[Model::CONFIG->total_time()];
  ys3 = new double[Model::CONFIG->total_time()];
  for (int i = 0; i < Model::CONFIG->total_time(); ++i) {
    xs1[i] = i;
    ys1[i] = i;
    ys2[i] = i;
    ys3[i] = i;
  }
}

void GPU::PlotLines::update_plot() {
  ImGui::Begin("Plots");
  t += ImGui::GetIO().DeltaTime;
  sdata1.AddPoint(t, Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()[6]);
  static int history = 10;
  if (ImPlot::BeginPlot("PfPR", ImVec2(-1,300))) {
    ImPlot::SetupAxes("Day","PfPR%",ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_RangeFit);
    ImPlot::SetupAxisLimits(ImAxis_X1,t - history, t, ImGuiCond_Always);
    ImPlot::PlotLine("Data 0", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), 0, sdata1.Offset, 2*sizeof(float));
    ImPlot::EndPlot();
  }
//  if(Model::GPU_SCHEDULER->current_time() % 30 == 0){
//    ys1[Model::GPU_SCHEDULER->current_time()] = Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()[6];
//    ys2[Model::GPU_SCHEDULER->current_time()] = Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()[7];
//    ys3[Model::GPU_SCHEDULER->current_time()] = Model::GPU_DATA_COLLECTOR->current_tf_by_therapy()[8];
//  }
//  if (ImPlot::BeginPlot("Line Plots", ImVec2(-1,300))) {
//    ImPlot::SetupAxes("Day","TF",ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_RangeFit);
//    ImPlot::SetupAxisLimits(ImAxis_X1,0,Model::CONFIG->total_time());
//    ImPlot::PlotLine("6", xs1, ys1, 3650);
//    ImPlot::PlotLine("7", xs1, ys2, 3650);
//    ImPlot::PlotLine("8", xs1, ys3, 3650);
//    ImPlot::EndPlot();
//  }
  ImGui::End();
}
