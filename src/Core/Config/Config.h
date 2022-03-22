/*
 * File:   Config.h
 * Author: nguyentran
 *
 * Created on March 27, 2013, 10:38 AM
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <date/date.h>

#include <string>
#include <vector>

#include "ConfigItem.h"
#include "Core/MultinomialDistributionGenerator.h"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "CustomConfigItem.h"
#include "Spatial/Location.h"

class Model;

class Config {
  DISALLOW_COPY_AND_ASSIGN(Config)

  DISALLOW_MOVE(Config)

public:
  POINTER_PROPERTY(Model, model)

  std::vector<IConfigItem *> config_items {};

  CONFIG_ITEM(starting_date, date::year_month_day, date::year_month_day { date::year { 1999 } / 1 / 1 })
  CONFIG_ITEM(ending_date, date::year_month_day, date::year_month_day { date::year { 1999 } / 1 / 2 })

  CUSTOM_CONFIG_ITEM(total_time, 100)

  CONFIG_ITEM(start_collect_data_day, int, 365)

  CONFIG_ITEM(number_of_tracking_days, int, 0)
  CONFIG_ITEM(p_infection_from_an_infectious_bite, double, 0.0)

  CONFIG_ITEM(age_structure, std::vector<int>, std::vector<int> {})
  CONFIG_ITEM(initial_age_structure, std::vector<int>, std::vector<int> {})

  CONFIG_ITEM(tf_testing_day, int, 28)

  CONFIG_ITEM(days_to_clinical_under_five, int, 4)
  CONFIG_ITEM(days_to_clinical_over_five, int, 6)
  CONFIG_ITEM(days_mature_gametocyte_under_five, int, 4)
  CONFIG_ITEM(days_mature_gametocyte_over_five, int, 6)

  CONFIG_ITEM(p_compliance, double, 1.0)
  CONFIG_ITEM(min_dosing_days, int, 1)

  CONFIG_ITEM(gametocyte_level_full, double, 1.0)

  CONFIG_ITEM(p_relapse, double, 0.01)
  CONFIG_ITEM(relapse_duration, int, 30)

  CONFIG_ITEM(allow_new_coinfection_to_cause_symtoms, bool, true)
  CONFIG_ITEM(update_frequency, int, 7)
  CONFIG_ITEM(report_frequency, int, 30)

  CONFIG_ITEM(tf_rate, double, 0.1)
  CONFIG_ITEM(using_free_recombination, bool, true)
  CONFIG_ITEM(tf_window_size, int, 60)

  CONFIG_ITEM(using_age_dependent_bitting_level, bool, false)
  CONFIG_ITEM(using_variable_probability_infectious_bites_cause_infection, bool, false)

  CONFIG_ITEM(fraction_mosquitoes_interrupted_feeding, double, 0.0)
  CONFIG_ITEM(inflation_factor, double, 0.01)

  CONFIG_ITEM(location_db, std::vector<Spatial::Location>,
              std::vector<Spatial::Location> { Spatial::Location(0, 0, 0, 10000) })

  CONFIG_ITEM(birth_rate, double, 0)

  CONFIG_ITEM(as_iov, double, 0.2)

  CONFIG_ITEM(death_rate_by_age_class, DoubleVector, DoubleVector())

  CONFIG_ITEM(mortality_when_treatment_fail_by_age_class, DoubleVector, DoubleVector())

  CONFIG_ITEM(parasite_density_level, ParasiteDensityLevel, ParasiteDensityLevel())

  CONFIG_ITEM(relative_infectivity, RelativeInfectivity, RelativeInfectivity())

  CONFIG_ITEM(pf_genotype_info, PfGeneInfo, PfGeneInfo())

  CONFIG_ITEM(initial_strategy_id, int, -1)

  CONFIG_ITEM(age_bracket_prob_individual_present_at_mda, IntVector, IntVector())

  CONFIG_ITEM(mean_prob_individual_present_at_mda, DoubleVector, DoubleVector())

  CONFIG_ITEM(sd_prob_individual_present_at_mda, DoubleVector, DoubleVector())

  CONFIG_ITEM(mda_therapy_id, int, 0)

  CONFIG_ITEM(artificial_rescaling_of_population_size, double, 1.0)

  CONFIG_ITEM(override_ec50_patterns, OverrideEC50Patterns, OverrideEC50Patterns())

  CONFIG_ITEM(mutation_mask, std::string, "")

  CONFIG_ITEM(mosquito_config, MosquitoConfig, MosquitoConfig())

  CONFIG_ITEM(within_chromosome_recombination_rate, double, 0)

  CUSTOM_CONFIG_ITEM(start_of_comparison_period, 0)

  CUSTOM_CONFIG_ITEM(number_of_age_classes, 0)

  CUSTOM_CONFIG_ITEM(number_of_locations, 1)

  CUSTOM_CONFIG_ITEM(spatial_distance_matrix, DoubleVector2())

  CUSTOM_CONFIG_ITEM(seasonal_info, SeasonalInfo())

  CUSTOM_CONFIG_ITEM(spatial_model, nullptr)

  CUSTOM_CONFIG_ITEM(immune_system_information, ImmuneSystemInformation())

  CUSTOM_CONFIG_ITEM(drug_db, nullptr)

  CUSTOM_CONFIG_ITEM(circulation_info, RelativeMovingInformation())

  CUSTOM_CONFIG_ITEM(relative_bitting_info, RelativeBittingInformation())

  CUSTOM_CONFIG_ITEM(therapy_db, TherapyPtrVector())

  CUSTOM_CONFIG_ITEM(strategy_db, StrategyPtrVector())

  CUSTOM_CONFIG_ITEM(initial_parasite_info, std::vector<InitialParasiteInfo>())

  CUSTOM_CONFIG_ITEM(preconfig_population_events, std::vector<Event *>())

  CUSTOM_CONFIG_ITEM(bitting_level_generator, MultinomialDistributionGenerator())
  CUSTOM_CONFIG_ITEM(moving_level_generator, MultinomialDistributionGenerator())

  CUSTOM_CONFIG_ITEM(prob_individual_present_at_mda_distribution, std::vector<beta_distribution_params>())

  VIRTUAL_PROPERTY_REF(double, modified_mutation_factor)

  VIRTUAL_PROPERTY_REF(double, modified_drug_half_life)

  VIRTUAL_PROPERTY_REF(double, modified_daily_cost_of_resistance)

  VIRTUAL_PROPERTY_REF(double, modified_mutation_probability)

public:
  GenotypeDatabase genotype_db {};

public:
  explicit Config(Model *model = nullptr);

  virtual ~Config();

  void read_from_file(const std::string &config_file_name = "config.yml");

  virtual void update_mutation_mask(const std::string &new_mask);
};

#endif /* CONFIG_H */
