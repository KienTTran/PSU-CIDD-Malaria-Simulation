//
// Created by Nguyen Tran on 3/17/2018.
//

#include <sstream>
#include "NestedMFTMultiLocationStrategy.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Population/Person.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Gpu/Therapies/Therapy.cuh"


GPU::NestedMFTMultiLocationStrategy::NestedMFTMultiLocationStrategy() : GPU::IStrategy(
    "NestedMFTMultiLocationStrategy", NestedMFTMultiLocation
) { }

GPU::NestedMFTMultiLocationStrategy::~NestedMFTMultiLocationStrategy() = default;

void GPU::NestedMFTMultiLocationStrategy::add_strategy(GPU::IStrategy* strategy) {
  strategy_list.push_back(strategy);
}

void GPU::NestedMFTMultiLocationStrategy::add_therapy(GPU::Therapy* therapy) { }

GPU::Therapy* GPU::NestedMFTMultiLocationStrategy::get_therapy(GPU::Person* person) {
  const auto loc = person->location();
  const auto p = Model::RANDOM->random_flat(0.0, 1.0);

  double sum = 0;
  for (auto i = 0; i < distribution[loc].size(); i++) {
    sum += distribution[loc][i];
    if (p <= sum) {
      return strategy_list[i]->get_therapy(person);
    }
  }
  return strategy_list[strategy_list.size() - 1]->get_therapy(person);
}

std::string GPU::NestedMFTMultiLocationStrategy::to_string() const {
  std::stringstream sstm;
  sstm << id << "-" << name << std::endl;

  for (auto i : distribution[0]) {
    sstm << i << ",";
  }
  sstm << std::endl;

  for (auto i : start_distribution[0]) {
    sstm << i << ",";
  }
  sstm << std::endl;
  return sstm.str();
}

void GPU::NestedMFTMultiLocationStrategy::update_end_of_time_step() {
  // update each strategy in the nest
  for (auto& strategy : strategy_list) {
    strategy->update_end_of_time_step();
  }
}

void GPU::NestedMFTMultiLocationStrategy::adjust_distribution(const int& time) {
  if (peak_after == -1) {
    // inflation every year
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
      const auto d_act = distribution[loc][0] * (1 + Model::CONFIG->inflation_factor() / 12);
      distribution[loc][0] = d_act;
      const auto other_d = (1 - d_act) / (distribution[loc].size() - 1);
      for (auto i = 1; i < distribution[loc].size(); i++) {
        distribution[loc][i] = other_d;
      }
    }
  } else {
    // increasing linearly
    if (time <= starting_time + peak_after) {
      if (distribution[0][0] < 1) {
        for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
          for (auto i = 0; i < distribution[loc].size(); i++) {

            auto dist = peak_after == 0 ? peak_distribution[loc][i] :
                        (peak_distribution[loc][i] - start_distribution[loc][i]) * (time - starting_time) /
                        peak_after + start_distribution[loc][i];
            dist = dist > peak_distribution[loc][i] ? peak_distribution[loc][i] : dist;
            distribution[loc][i] = dist;
          }
        }
      }
    }
  } //    std::cout << to_string() << std::endl;
}

void GPU::NestedMFTMultiLocationStrategy::adjust_started_time_point(const int& current_time) {
  starting_time = current_time;
  // update each strategy in the nest
  for (auto* strategy : strategy_list) {
    strategy->adjust_started_time_point(current_time);
  }
}

void GPU::NestedMFTMultiLocationStrategy::monthly_update() {
  adjust_distribution(Model::GPU_SCHEDULER->current_time());

  for (auto* strategy : strategy_list) {
    strategy->monthly_update();
  }

  // for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
  //   std::cout << distribution[loc] << std::endl;
  // }

}
