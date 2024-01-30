/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StrategyBuilder.cpp
 * Author: Merlin
 * 
 * Created on August 23, 2017, 11:03 AM
 */

#include "StrategyBuilder.cuh"
#include "IStrategy.cuh"
#include "SFTStrategy.cuh"
#include "Core/Config/Config.h"
#include "CyclingStrategy.cuh"
#include "AdaptiveCyclingStrategy.cuh"
#include "MFTStrategy.cuh"
#include "NestedMFTStrategy.cuh"
#include "MFTRebalancingStrategy.cuh"
#include "MFTMultiLocationStrategy.cuh"
#include "NestedMFTMultiLocationStrategy.cuh"
#include "NovelDrugIntroductionStrategy.cuh"

GPU::StrategyBuilder::StrategyBuilder() = default;

GPU::StrategyBuilder::~StrategyBuilder() = default;

GPU::IStrategy* GPU::StrategyBuilder::build(const YAML::Node &ns, const int &strategy_id, Config* config) {
  const auto type = GPU::IStrategy::StrategyTypeMap[ns["type"].as<std::string>()];
  switch (type) {
    case GPU::IStrategy::SFT:
      return buildSFTStrategy(ns, strategy_id, config);
    case GPU::IStrategy::Cycling:
      return buildCyclingStrategy(ns, strategy_id, config);
    case GPU::IStrategy::AdaptiveCycling:
      return buildAdaptiveCyclingStrategy(ns, strategy_id, config);
    case GPU::IStrategy::MFT:
      return buildMFTStrategy(ns, strategy_id, config);
    case GPU::IStrategy::MFTRebalancing:
      return buildMFTRebalancingStrategy(ns, strategy_id, config);
    case GPU::IStrategy::NestedMFT:
      return buildNestedSwitchingStrategy(ns, strategy_id, config);
    case GPU::IStrategy::MFTMultiLocation:
      return buildMFTMultiLocationStrategy(ns, strategy_id, config);
    case GPU::IStrategy::NestedMFTMultiLocation:
      return buildNestedMFTDifferentDistributionByLocationStrategy(ns,
                                                                   strategy_id,
                                                                   config);
    case GPU::IStrategy::NovelDrugIntroduction:
      return buildNovelDrugIntroductionStrategy(ns, strategy_id, config);
    default:
      return nullptr;
  }
}

void GPU::StrategyBuilder::add_therapies(const YAML::Node &ns, GPU::IStrategy* result, Config* config) {
  for (auto i = 0; i < ns["therapy_ids"].size(); i++) {
    result->add_therapy(config->gpu_therapy_db()[ns["therapy_ids"][i].as<int>()]);
  }
}

void GPU::StrategyBuilder::add_distributions(const YAML::Node &ns, DoubleVector &v) {
  for (auto i = 0; i < ns.size(); i++) {
    v.push_back(ns[i].as<double>());
  }
}

GPU::IStrategy* GPU::StrategyBuilder::buildSFTStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::SFTStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();
  result->add_therapy(config->gpu_therapy_db()[ns["therapy_id"].as<int>()]);
  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildCyclingStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::CyclingStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->cycling_time = ns["cycling_time"].as<int>();
  result->next_switching_day = ns["cycling_time"].as<int>();

  add_therapies(ns, result, config);

  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildAdaptiveCyclingStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::AdaptiveCyclingStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->trigger_value = ns["trigger_value"].as<double>();
  result->delay_until_actual_trigger = ns["delay_until_actual_trigger"].as<int>();
  result->turn_off_days = ns["turn_off_days"].as<int>();

  add_therapies(ns, result, config);
  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildMFTStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::MFTStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["distribution"], result->distribution);
  add_therapies(ns, result, config);
  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildNestedSwitchingStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::NestedMFTStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["start_distribution"], result->start_distribution);
  add_distributions(ns["start_distribution"], result->distribution);
  add_distributions(ns["peak_distribution"], result->peak_distribution);

  result->peak_after = ns["peak_after"].as<int>();

  for (int i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(
        config->gpu_strategy_db()[ns["strategy_ids"][i].as<int>()]);
  }

  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildMFTRebalancingStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::MFTRebalancingStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["distribution"], result->distribution);
  add_distributions(ns["distribution"], result->next_distribution);

  add_therapies(ns, result, config);

  result->update_duration_after_rebalancing = ns["update_duration_after_rebalancing"].as<int>();
  result->delay_until_actual_trigger = ns["delay_until_actual_trigger"].as<int>();
  result->latest_adjust_distribution_time = 0;

  return result;
}

GPU::IStrategy*
GPU::StrategyBuilder::buildMFTMultiLocationStrategy(const YAML::Node &ns, const int &strategy_id, Config* config) {
  auto* result = new GPU::MFTMultiLocationStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->distribution.clear();
  result->distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  result->start_distribution.clear();
  result->start_distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  result->peak_distribution.clear();
  result->peak_distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["start_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["start_distribution"][input_loc], result->distribution[loc]);
  }
  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["start_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["start_distribution"][input_loc], result->start_distribution[loc]);
  }

  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["peak_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["peak_distribution"][input_loc], result->peak_distribution[loc]);
  }

  add_therapies(ns, result, config);
  result->peak_after = ns["peak_after"].as<int>();
  return result;
}

GPU::IStrategy* GPU::StrategyBuilder::buildNestedMFTDifferentDistributionByLocationStrategy(const YAML::Node &ns,
                                                                                  const int &strategy_id,
                                                                                  Config* config) {
  auto* result = new GPU::NestedMFTMultiLocationStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->distribution.clear();
  result->distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  result->start_distribution.clear();
  result->start_distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  result->peak_distribution.clear();
  result->peak_distribution.resize(static_cast<unsigned long long int>(config->number_of_locations()));

  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["start_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["start_distribution"][input_loc], result->distribution[loc]);
  }
  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["start_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["start_distribution"][input_loc], result->start_distribution[loc]);
  }

  for (auto loc = 0; loc < config->number_of_locations(); loc++) {
    auto input_loc = ns["peak_distribution"].size() < config->number_of_locations() ? 0 : loc;
    add_distributions(ns["peak_distribution"][input_loc], result->peak_distribution[loc]);
  }

  for (auto i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(config->gpu_strategy_db()[ns["strategy_ids"][i].as<int>()]);
  }

  result->peak_after = ns["peak_after"].as<int>();
  //    std::cout << result->to_string() << std::endl;

  return result;
}

GPU::IStrategy*
GPU::StrategyBuilder::buildNovelDrugIntroductionStrategy(const YAML::Node &ns, const int strategy_id, Config* config) {
  auto* result = new GPU::NovelDrugIntroductionStrategy();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["start_distribution"], result->start_distribution);
  add_distributions(ns["start_distribution"], result->distribution);
  add_distributions(ns["peak_distribution"], result->peak_distribution);

  result->peak_after = ns["peak_after"].as<int>();

  for (int i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(
        config->gpu_strategy_db()[ns["strategy_ids"][i].as<int>()]);
  }

  result->newly_introduced_strategy_id = ns["newly_introduced_strategy_id"].as<int>();
  result->tf_threshold = ns["tf_threshold"].as<double>();

  result->replacement_fraction =  ns["replacement_fraction"].as<double>();
  result->replacement_duration =  ns["replacement_duration"].as<double>();

  return result;
}
