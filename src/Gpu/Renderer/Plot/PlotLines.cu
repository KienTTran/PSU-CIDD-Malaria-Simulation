//
// Created by ktt on 2/20/24.
//
#include "PlotLines.cuh"
#include <imgui_stdlib.h>
#include "Model.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/MDC/ModelDataCollector.cuh"

GPU::PlotLines::~PlotLines() {

}

GPU::PlotLines::PlotLines() {

}

void GPU::PlotLines::start_plot() {
}

void GPU::PlotLines::update_plot() {
  ImGui::BulletText("This example assumes 60 FPS. Higher FPS requires larger buffer size.");
  static ScrollingBuffer sdata1, sdata2;
  ImVec2 mouse = ImGui::GetMousePos();
  static float t = 0;
  t += Model::GPU_SCHEDULER->current_time()*1.0f;
  sdata1.AddPoint(t, Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[0] * 100.0f);
  sdata2.AddPoint(t, Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[1] * 100.0f);

  static int history = 3650;
//  ImGui::SliderFloat("History",&history,1,30,"%.1f s");

  static ImPlotAxisFlags flags = ImPlotAxisFlags_NoTickLabels;
  if (ImPlot::BeginPlot("##Scrolling", ImVec2(-1,150))) {
    ImPlot::SetupAxes(nullptr, nullptr, flags, flags);
    ImPlot::SetupAxisLimits(ImAxis_X1,t - history, t, ImGuiCond_Always);
    ImPlot::SetupAxisLimits(ImAxis_Y1,0,100);
    ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL,0.5f);
    ImPlot::PlotLine("Pop 0", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), 0, sdata1.Offset, 2 * sizeof(int));
    ImPlot::PlotLine("Pop 1", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), 0, sdata2.Offset, 2*sizeof(int));
    ImPlot::EndPlot();
  }
}