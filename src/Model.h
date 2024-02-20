/* 
 * File:   Model.h
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:26 PM
 */

#ifndef MODEL_H
#define    MODEL_H

#include <vector>
#include "Core/Scheduler.h"
#include "Gpu/Core/Scheduler.cuh"
#include "Core/PropertyMacro.h"
#include "Malaria/ITreatmentCoverageModel.h"
#include "Population/ClinicalUpdateFunction.h"
#include "Population/ImmunityClearanceUpdateFunction.h"
#include "Gpu/Population/ClinicalUpdateFunction.cuh"
#include "Gpu/Population/ImmunityClearanceUpdateFunction.cuh"

namespace GPU{
    class Random;
    class PopulationKernel;
    class Renderer;
    class RenderEntity;
    class Utils;
    class Mosquito;
    class Population;
    class IStrategy;
    class ModelDataCollector;
    class Reporter;
    class Scheduler;
    class Plot;
}

class Mosquito;

class Scheduler;

class Population;

class Config;

class ModelDataCollector;

class Reporter;

class Random;

class Model {
 DISALLOW_COPY_AND_ASSIGN(Model)
 void end_time_step();

 DISALLOW_MOVE(Model)

 POINTER_PROPERTY(Config, config)

 POINTER_PROPERTY(::Scheduler, scheduler)

 POINTER_PROPERTY(GPU::Scheduler, gpu_scheduler)

 POINTER_PROPERTY(::Population, population)

 POINTER_PROPERTY(ClinicalUpdateFunction, progress_to_clinical_update_function)
 POINTER_PROPERTY(ImmunityClearanceUpdateFunction, immunity_clearance_update_function)
 POINTER_PROPERTY(ImmunityClearanceUpdateFunction, having_drug_update_function)
 POINTER_PROPERTY(ImmunityClearanceUpdateFunction, clinical_update_function)

 POINTER_PROPERTY(::Random, random)

 POINTER_PROPERTY(ModelDataCollector, data_collector)

 PROPERTY_REF(std::vector<GPU::Reporter *>, reporters)

 PROPERTY_REF(unsigned long, initial_seed_number)

 PROPERTY_REF(std::string, config_filename)

 PROPERTY_REF(int, cluster_job_number)

 PROPERTY_REF(std::string, output_path)

 PROPERTY_REF(std::string, tme_filename)

 PROPERTY_REF(std::string, override_parameter_filename)

 PROPERTY_REF(int, override_parameter_line_number) // base 1
 PROPERTY_REF(int, gui_type)

 PROPERTY_REF(bool, is_farm_output)

 PROPERTY_REF(std::string, reporter_type)

POINTER_PROPERTY(::Mosquito, mosquito)

POINTER_PROPERTY(GPU::ModelDataCollector, gpu_data_collector)
POINTER_PROPERTY(GPU::Population, gpu_population)
POINTER_PROPERTY(GPU::Renderer, gpu_renderer)
POINTER_PROPERTY(GPU::RenderEntity, gpu_render_entity)
POINTER_PROPERTY(GPU::Utils, gpu_utils)
POINTER_PROPERTY(GPU::Random, gpu_random)
POINTER_PROPERTY(GPU::PopulationKernel, gpu_population_kernel)
POINTER_PROPERTY(GPU::Mosquito, gpu_mosquito)
POINTER_PROPERTY(GPU::ClinicalUpdateFunction, gpu_progress_to_clinical_update_function)
POINTER_PROPERTY(GPU::ImmunityClearanceUpdateFunction, gpu_immunity_clearance_update_function)
POINTER_PROPERTY(GPU::ImmunityClearanceUpdateFunction, gpu_having_drug_update_function)
POINTER_PROPERTY(GPU::ImmunityClearanceUpdateFunction, gpu_clinical_update_function)
POINTER_PROPERTY(GPU::Plot, gpu_plot)
 public:
  static Model *MODEL;
  static Config *CONFIG;
  static ::Random *RANDOM;
  static ::Scheduler *SCHEDULER;
  static ::Population *POPULATION;
  static ::ModelDataCollector *DATA_COLLECTOR;
  static ::Mosquito *MOSQUITO;

  static GPU::Scheduler *GPU_SCHEDULER;
  static GPU::ModelDataCollector *GPU_DATA_COLLECTOR;
  static GPU::Population *GPU_POPULATION;
  static GPU::PopulationKernel *GPU_POPULATION_KERNEL;
  static GPU::Renderer *GPU_RENDERER;
  static GPU::RenderEntity *GPU_RENDER_ENTITY;
  static GPU::Utils *GPU_UTILS;
  static GPU::Random *GPU_RANDOM;
  static GPU::Mosquito *GPU_MOSQUITO;
  static GPU::Plot *GPU_PLOT;
  bool model_finished = false;

  static IStrategy *TREATMENT_STRATEGY;
  static GPU::IStrategy *GPU_TREATMENT_STRATEGY;
  static ITreatmentCoverageModel *TREATMENT_COVERAGE;
  // static std::shared_ptr<spdlog::logger> LOGGER;

  explicit Model(const int &object_pool_size = 100000);

  virtual ~Model();

  void set_treatment_strategy(const int &strategy_id);

  void set_treatment_coverage(ITreatmentCoverageModel *tcm);

  void build_initial_treatment_coverage();

  void initialize();

  static void initialize_object_pool(const int &size = 100000);

  static void release_object_pool();

  void before_run();

  void run();

  void after_run();

  void begin_time_step();

  void daily_update();

  void monthly_update();

  void yearly_update();

  void release();

  void report_begin_of_time_step();

  void monthly_report();

  void add_reporter(GPU::Reporter *reporter);

  double get_seasonal_factor(const date::sys_days &today, const int &location) const;

 private:
    IStrategy *treatment_strategy_{nullptr};
    GPU::IStrategy *gpu_treatment_strategy_{nullptr};
    ITreatmentCoverageModel *treatment_coverage_{nullptr};
};

#endif    /* MODEL_H */
