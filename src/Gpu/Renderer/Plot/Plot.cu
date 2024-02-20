//
// Created by ktt on 2/20/24.
//

#include "Plot.cuh"

GPU::Plot::Plot() {
  plots_ = TVector<GPU::Plot*>();
}

GPU::Plot::~Plot() {
  for(auto plt : plots_) {
    delete plt;
  }
  plots_.clear();
}

void GPU::Plot::init() {
}

void GPU::Plot::add_plot(GPU::Plot* plot) {
  plots_.push_back(plot);
}

void GPU::Plot::start_plot() {
  for(auto plt : plots_) {
    plt->start_plot();
  }
}

void GPU::Plot::update_plot() {
  for(auto plt : plots_) {
    plt->update_plot();
  }
}

