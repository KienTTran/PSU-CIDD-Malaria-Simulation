#include "NestedMFTStrategy.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Therapies/Therapy.cuh"

void GPU::NestedMFTStrategy::add_strategy(GPU::IStrategy* strategy) {
  strategy_list.push_back(strategy);
}

void GPU::NestedMFTStrategy::add_therapy(GPU::Therapy* therapy) { }

GPU::Therapy* GPU::NestedMFTStrategy::get_therapy(GPU::Person* person) {
  const auto p = Model::RANDOM->random_flat(0.0, 1.0);

  double sum = 0;
  for (auto i = 0; i < distribution.size(); i++) {
    sum += distribution[i];
    if (p <= sum) {
      return strategy_list[i]->get_therapy(person);
    }
  }
  return strategy_list[strategy_list.size() - 1]->get_therapy(person);
}

std::string GPU::NestedMFTStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << "-";

  for (auto i : distribution) {
    sstm << i << "::";
  }
  return sstm.str();
}

void GPU::NestedMFTStrategy::adjust_started_time_point(const int& current_time) {
  // update each strategy in the nest
  for (auto& strategy : strategy_list) {
    strategy->adjust_started_time_point(current_time);
  }
  starting_time = current_time;
}

void GPU::NestedMFTStrategy::update_end_of_time_step() {
  // update each strategy in the nest
  for (auto* strategy : strategy_list) {
    strategy->update_end_of_time_step();
  }
}

void GPU::NestedMFTStrategy::monthly_update() {
  adjust_distribution(Model::GPU_SCHEDULER->current_time());

  for (auto* strategy : strategy_list) {
    strategy->monthly_update();
  }
//  std::cout << distribution[0] << "-" << distribution[1] << std::endl;
}

void GPU::NestedMFTStrategy::adjust_distribution(const int& time) {
  if (time <= starting_time + peak_after) {
    for (auto i = 0; i < distribution.size(); i++) {
      const auto dist = (peak_distribution[i] - start_distribution[i]) * (time - starting_time) / peak_after +
                        start_distribution[i];
      distribution[i] = dist;
    }
  } else {
    if (peak_distribution[0] != distribution[0]) {
      for (auto i = 0; i < distribution.size(); i++) {
        distribution[i] = peak_distribution[i];
      }
    }
  }
}
