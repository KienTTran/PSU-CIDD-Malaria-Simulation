/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StrategyBuilder.h
 * Author: Merlin
 *
 * Created on August 23, 2017, 11:03 AM
 */

#ifndef STRATEGYBUILDER_CUH
#define STRATEGYBUILDER_CUH

#include <yaml-cpp/yaml.h>

#include "Core/TypeDef.h"
#include "Core/PropertyMacro.h"

namespace GPU{
    class IStrategy;
    class StrategyBuilder;
}

class Config;

class GPU::StrategyBuilder {
 DISALLOW_COPY_AND_ASSIGN(StrategyBuilder)

 public:
  StrategyBuilder();

  virtual ~StrategyBuilder();

  static GPU::IStrategy *build(const YAML::Node &ns, const int &strategy_id, Config *config);

  static void add_therapies(const YAML::Node &ns, GPU::IStrategy *result, Config *config);

  static void add_distributions(const YAML::Node &ns, DoubleVector &v);

  static GPU::IStrategy *buildSFTStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildCyclingStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildAdaptiveCyclingStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildMFTStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildMFTRebalancingStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildNestedSwitchingStrategy(const YAML::Node &ns, const int &strategy_id, Config *config);

  static GPU::IStrategy *buildMFTMultiLocationStrategy(const YAML::Node &node, const int &id, Config *config);

  static GPU::IStrategy *
  buildNestedMFTDifferentDistributionByLocationStrategy(const YAML::Node &ns, const int &strategy_id,
                                                        Config *config);

  static GPU::IStrategy *buildNovelDrugIntroductionStrategy(const YAML::Node &ns, const int strategy_id, Config *config);
};

#endif /* STRATEGYBUILDER_CUH */
