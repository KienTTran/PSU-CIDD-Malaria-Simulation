/*
 * File:   Model.cpp
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:26 PM
 */
#include "Model.h"
#include <fmt/format.h>

#include "Constants.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Events/BirthdayEvent.h"
#include "Events/CirculateToTargetLocationNextDayEvent.h"
#include "Events/EndClinicalByNoTreatmentEvent.h"
#include "Events/EndClinicalDueToDrugResistanceEvent.h"
#include "Events/EndClinicalEvent.h"
#include "Events/MatureGametocyteEvent.h"
#include "Events/MoveParasiteToBloodEvent.h"
#include "Events/Population/ImportationEvent.h"
#include "Events/Population/ImportationPeriodicallyEvent.h"
#include "Events/ProgressToClinicalEvent.h"
#include "Events/ReturnToResidenceEvent.h"
#include "Events/SwitchImmuneComponentEvent.h"
#include "Events/TestTreatmentFailureEvent.h"
#include "Events/UpdateWhenDrugIsPresentEvent.h"
#include "Helpers/ObjectHelpers.h"
#include "Helpers/TimeHelpers.h"
#include "MDC/ModelDataCollector.h"
#include "Malaria/SteadyTCM.h"
#include "Mosquito/Mosquito.h"
#include "Population/ClonalParasitePopulation.h"
#include "Population/ImmuneSystem.h"
#include "Population/Person.h"
#include "Population/Population.h"
#include "Population/SingleHostClonalParasitePopulations.h"
#include "Reporters/Reporter.h"
#include "Strategies/IStrategy.h"
#include "Therapies/Drug.h"
#include "easylogging++.h"
#include "Spatial/SpatialModel.hxx"
#include "Gpu/Renderer/Renderer.h"
#include "Gpu/RenderEntity.cuh"
#include "Gpu/Utils.cuh"
#include "Gpu/Random.cuh"
#include "Gpu/Population/Population.cuh"

Model* Model::MODEL = nullptr;
Config* Model::CONFIG = nullptr;
Random* Model::RANDOM = nullptr;
Scheduler* Model::SCHEDULER = nullptr;
ModelDataCollector* Model::DATA_COLLECTOR = nullptr;
Population* Model::POPULATION = nullptr;
Mosquito* Model::MOSQUITO = nullptr;
IStrategy* Model::TREATMENT_STRATEGY = nullptr;
ITreatmentCoverageModel* Model::TREATMENT_COVERAGE = nullptr;
// std::shared_ptr<spdlog::logger> LOGGER;
Renderer* Model::RENDERER = nullptr;
GPU::RenderEntity* Model::GPU_RENDER_ENTITY = nullptr;
GPU::Random* Model::GPU_RANDOM = nullptr;
GPU::Utils* Model::GPU_UTILS = nullptr;
GPU::Population* Model::GPU_POPULATION = nullptr;

Model::Model(const int& object_pool_size) {
  initialize_object_pool(object_pool_size);
  random_ = new Random();
  config_ = new Config(this);
  scheduler_ = new Scheduler(this);
  population_ = new Population(this);
  //  mosquito = new Mosquito(this);
  data_collector_ = new ModelDataCollector(this);
  mosquito = new Mosquito(this);
  renderer_ = new Renderer(this);
  gpu_render_entity_ = new GPU::RenderEntity(this);
  gpu_random_ = new GPU::Random();
  gpu_utils_ = new GPU::Utils();
  gpu_population_ = new GPU::Population();

  MODEL = this;
  CONFIG = config_;
  SCHEDULER = scheduler_;
  RANDOM = random_;
  DATA_COLLECTOR = data_collector_;
  POPULATION = population_;
  MOSQUITO = mosquito;

  RENDERER = renderer_;
  GPU_RENDER_ENTITY = gpu_render_entity_;
  GPU_RANDOM = gpu_random_;
  GPU_UTILS = gpu_utils_;
  GPU_POPULATION = gpu_population_;

  // LOGGER = spdlog::stdout_logger_mt("console");

  progress_to_clinical_update_function_ = new ClinicalUpdateFunction(this);
  immunity_clearance_update_function_ = new ImmunityClearanceUpdateFunction(this);
  having_drug_update_function_ = new ImmunityClearanceUpdateFunction(this);
  clinical_update_function_ = new ImmunityClearanceUpdateFunction(this);

  reporters_ = std::vector<Reporter*>();

  initial_seed_number_ = 0;
  config_filename_ = "config.yml";
  tme_filename_ = "tme.txt";
  override_parameter_filename_ = "";
  override_parameter_line_number_ = -1;
  gui_type_ = -1;
  is_farm_output_ = false;
  cluster_job_number_ = 0;
  reporter_type_ = "";
}

Model::~Model() {
  release();

  release_object_pool();
}

void Model::set_treatment_strategy(const int& strategy_id) {
  treatment_strategy_ = strategy_id == -1 ? nullptr : config_->strategy_db()[strategy_id];
  TREATMENT_STRATEGY = treatment_strategy_;

  treatment_strategy_->adjust_started_time_point(Model::SCHEDULER->current_time());

  //
  // if (treatment_strategy_->get_type() == IStrategy::NestedSwitching) {
  //   dynamic_cast<NestedSwitchingStrategy *>(treatment_strategy_)->initialize_update_time(config_);
  // }
  // if (treatment_strategy_->get_type() == IStrategy::NestedSwitchingDifferentDistributionByLocation) {
  //   dynamic_cast<NestedMFTMultiLocationStrategy *>(treatment_strategy_)->initialize_update_time(config_);
  // }
}

void Model::set_treatment_coverage(ITreatmentCoverageModel* tcm) {
  if (treatment_coverage_ != tcm) {
    if (tcm->p_treatment_less_than_5.empty() || tcm->p_treatment_more_than_5.empty()) {
      // copy current value
      tcm->p_treatment_less_than_5 = treatment_coverage_->p_treatment_less_than_5;
      tcm->p_treatment_more_than_5 = treatment_coverage_->p_treatment_more_than_5;
    }

    ObjectHelpers::delete_pointer<ITreatmentCoverageModel>(treatment_coverage_);
  }
  treatment_coverage_ = tcm;
  TREATMENT_COVERAGE = tcm;
}

void Model::build_initial_treatment_coverage() {
  auto* tcm = new SteadyTCM();
  for (auto& location : config_->location_db()) {
    tcm->p_treatment_less_than_5.push_back(location.p_treatment_less_than_5);
    tcm->p_treatment_more_than_5.push_back(location.p_treatment_more_than_5);
  }
  set_treatment_coverage(tcm);
}

void Model::initialize() {
  LOG(INFO) << "Model initilizing...";

  LOG(INFO) << fmt::format("Read input file: {}", config_filename_);
  // Read input file
  config_->read_from_file(config_filename_);

  auto now = std::chrono::high_resolution_clock::now();
  auto milliseconds = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());

  LOG(INFO) << "Initialize Random";
  // Initialize Random Seed
  initial_seed_number_ = Model::CONFIG->initial_seed_number() <= 0 ? static_cast<unsigned long>(milliseconds.count()) : Model::CONFIG->initial_seed_number();

  LOG(INFO) << "Initialize Random";
  // Initialize Random Seed
  random_->initialize(initial_seed_number_);

  // add reporter here
  if (reporter_type_.empty()) {
    add_reporter(Reporter::MakeReport(Reporter::MONTHLY_REPORTER));
  } else {
    if (Reporter::ReportTypeMap.find(reporter_type_) != Reporter::ReportTypeMap.end()) {
      add_reporter(Reporter::MakeReport(Reporter::ReportTypeMap[reporter_type_]));
    }
  }

  LOG(INFO) << "Initialing reports";
  // initialize reporters
  for (auto* reporter : reporters_) {
    reporter->initialize();
  }

  LOG(INFO) << "Initialzing scheduler";
  LOG(INFO) << "Starting day is " << CONFIG->starting_date();
  // initialize scheduler
  scheduler_->initialize(CONFIG->starting_date(), config_->total_time());

  LOG(INFO) << "Initialing initial strategy";
  // set treatment strategy
  set_treatment_strategy(config_->initial_strategy_id());

  LOG(INFO) << "Initialing initial treatment coverage model";
  build_initial_treatment_coverage();

  LOG(INFO) << "Initializing data collector";
  // initialize data_collector
  data_collector_->initialize();

  LOG(INFO) << "Initializing population";
  // initialize Population
  population_->initialize();

  LOG(INFO) << "Initializing movement model...";
  config_->spatial_model()->prepare();

  LOG(INFO) << "Initializing mosquito";
  // initialize Population
  mosquito->initialize(config_);

  LOG(INFO) << "Introducing initial cases";
  // initialize infected_cases
  population_->introduce_initial_cases();

  LOG(INFO) << "Initializing GPU::Population";
  gpu_population_->init();

  LOG(INFO) << "Initializing GPU::Random";
  //Random always init with max population size
  gpu_random_->init(Model::CONFIG->n_people_init()*Model::CONFIG->gpu_config().pre_allocated_mem_ratio,initial_seed_number_);

  LOG(INFO) << "Initializing GPU::Utils";
  gpu_utils_->init();

  LOG(INFO) << "Initializing GPU::RenderEntity";
  gpu_render_entity_->init_entity();//send h_population to render

  LOG(INFO) << "Initializing GPU::Renderer";
  renderer_->init(gpu_render_entity_);

  LOG(INFO) << "Schedule for population event (if configured)";
  for (auto* event : config_->preconfig_population_events()) {
    scheduler_->schedule_population_event(event);
  }
}

void Model::initialize_object_pool(const int& size) {
  BirthdayEvent::InitializeObjectPool(size);
  ProgressToClinicalEvent::InitializeObjectPool(size);
  EndClinicalDueToDrugResistanceEvent::InitializeObjectPool(size);
  UpdateWhenDrugIsPresentEvent::InitializeObjectPool(size);
  EndClinicalEvent::InitializeObjectPool(size);
  EndClinicalByNoTreatmentEvent::InitializeObjectPool(size);
  MatureGametocyteEvent::InitializeObjectPool(size);
  MoveParasiteToBloodEvent::InitializeObjectPool(size);
  CirculateToTargetLocationNextDayEvent::InitializeObjectPool(size);
  ReturnToResidenceEvent::InitializeObjectPool(size);
  SwitchImmuneComponentEvent::InitializeObjectPool(size);
  ImportationPeriodicallyEvent::InitializeObjectPool(size);
  ImportationEvent::InitializeObjectPool(size);
  TestTreatmentFailureEvent::InitializeObjectPool(size);

  ClonalParasitePopulation::InitializeObjectPool(size);
  SingleHostClonalParasitePopulations::InitializeObjectPool();

  Drug::InitializeObjectPool(size);
  DrugsInBlood::InitializeObjectPool(size);

  //    InfantImmuneComponent::InitializeObjectPool(size);
  //    NonInfantImmuneComponent::InitializeObjectPool(size);

  ImmuneSystem::InitializeObjectPool(size);
  Person::InitializeObjectPool(size);
}

void Model::release_object_pool() {
  //    std::cout << "Release object pool" << std::endl;
  Person::ReleaseObjectPool();
  ImmuneSystem::ReleaseObjectPool();

  // TODO: Investigate why?
  //    InfantImmuneComponent::ReleaseObjectPool();
  //    NonInfantImmuneComponent::ReleaseObjectPool();

  DrugsInBlood::ReleaseObjectPool();
  Drug::ReleaseObjectPool();

  SingleHostClonalParasitePopulations::ReleaseObjectPool();
  ClonalParasitePopulation::ReleaseObjectPool();

  TestTreatmentFailureEvent::ReleaseObjectPool();
  ImportationEvent::ReleaseObjectPool();
  ImportationPeriodicallyEvent::ReleaseObjectPool();
  SwitchImmuneComponentEvent::ReleaseObjectPool();
  ReturnToResidenceEvent::ReleaseObjectPool();
  CirculateToTargetLocationNextDayEvent::ReleaseObjectPool();
  MoveParasiteToBloodEvent::ReleaseObjectPool();
  MatureGametocyteEvent::ReleaseObjectPool();
  EndClinicalByNoTreatmentEvent::ReleaseObjectPool();
  EndClinicalEvent::ReleaseObjectPool();
  UpdateWhenDrugIsPresentEvent::ReleaseObjectPool();
  EndClinicalDueToDrugResistanceEvent::ReleaseObjectPool();
  ProgressToClinicalEvent::ReleaseObjectPool();
  BirthdayEvent::ReleaseObjectPool();
}

void Model::run() {
  LOG(INFO) << "Model starting...";
  before_run();
  if(Model::CONFIG->render_config().display_gui && Model::CONFIG->gpu_config().pre_allocated_mem_ratio > 1.0){
    std::thread scheduler_thread(&Scheduler::run, scheduler_);
    renderer_->start();
    scheduler_thread.join();
  }
  else{
    Model::CONFIG->render_config().display_gui = false;
    scheduler_->run();
  }
  after_run();
  LOG(INFO) << "Model finished.";
}

void Model::before_run() {
  LOG(INFO) << "Perform before run events";
  for (auto* reporter : reporters_) {
    reporter->before_run();
  }
}

void Model::after_run() {
  LOG(INFO) << "Perform after run events";

  data_collector_->update_after_run();

  for (auto* reporter : reporters_) {
    reporter->after_run();
  }
}

void Model::begin_time_step() {
  // reset daily variables
  data_collector_->begin_time_step();
  report_begin_of_time_step();
}

void Model::daily_update() {
  population_->update_all_individuals();
  // for safety remove all dead by calling perform_death_event
  population_->perform_death_event();
  population_->perform_birth_event();

  // update current foi should be call after perform death, birth event
  // in order to obtain the right all alive individuals,
  // infection event will use pre-calculated individual relative biting rate to infect new infections
  // circulation event will use pre-calculated individual relative moving rate to migrate individual to new location
  population_->update_current_foi();

  population_->perform_infection_event();
//  population_->perform_circulation_event();
  gpu_population_->perform_circulation_event();

  // infect new mosquito cohort in prmc must be run after population perform infection event and update current foi
  // because the prmc at the tracking index will be overridden with new cohort to use N days later and
  // infection event used the prmc at the tracking index for the today infection
  auto tracking_index = scheduler_->current_time() % config_->number_of_tracking_days();
  mosquito->infect_new_cohort_in_PRMC(config_, random_, population_, tracking_index);

  // this function must be called after mosquito infect new cohort in prmc
  population_->persist_current_force_of_infection_to_use_N_days_later();
}

void Model::monthly_update() {
  monthly_report();

  // reset monthly variables
  data_collector()->monthly_update();

  treatment_strategy_->monthly_update();

  treatment_coverage_->monthly_update();
}

// ReSharper disable once CppMemberFunctionMayBeConst
void Model::yearly_update() {
  data_collector_->perform_yearly_update();
}

void Model::end_time_step() {
  // update / calculate daily UTL
  data_collector_->end_of_time_step();

  // check to switch strategy
  treatment_strategy_->update_end_of_time_step();
}

void Model::release() {
  //    std::cout << "Model Release" << std::endl;
  ObjectHelpers::delete_pointer<ClinicalUpdateFunction>(progress_to_clinical_update_function_);
  ObjectHelpers::delete_pointer<ImmunityClearanceUpdateFunction>(immunity_clearance_update_function_);
  ObjectHelpers::delete_pointer<ImmunityClearanceUpdateFunction>(having_drug_update_function_);
  ObjectHelpers::delete_pointer<ImmunityClearanceUpdateFunction>(clinical_update_function_);

  ObjectHelpers::delete_pointer<Population>(population_);
  //   ObjectHelpers::DeletePointer<ExternalPopulation>(external_population_);
  ObjectHelpers::delete_pointer<Scheduler>(scheduler_);
  ObjectHelpers::delete_pointer<ModelDataCollector>(data_collector_);

  treatment_strategy_ = nullptr;
  ObjectHelpers::delete_pointer<ITreatmentCoverageModel>(treatment_coverage_);

  ObjectHelpers::delete_pointer<Config>(config_);
  ObjectHelpers::delete_pointer<Random>(random_);

  for (Reporter* reporter : reporters_) {
    ObjectHelpers::delete_pointer<Reporter>(reporter);
  }
  reporters_.clear();

  MODEL = nullptr;
  CONFIG = nullptr;
  SCHEDULER = nullptr;
  RANDOM = nullptr;
  DATA_COLLECTOR = nullptr;
  POPULATION = nullptr;
  TREATMENT_STRATEGY = nullptr;
  TREATMENT_COVERAGE = nullptr;
}

void Model::monthly_report() {
  data_collector_->perform_population_statistic();

  for (auto* reporter : reporters_) {
    reporter->monthly_report();
  }
}

void Model::report_begin_of_time_step() {
  for (auto* reporter : reporters_) {
    reporter->begin_time_step();
  }
}

void Model::add_reporter(Reporter* reporter) {
  reporters_.push_back(reporter);
  reporter->set_model(this);
}

//
//double Model::get_seasonal_factor(const date::sys_days& today, const int& location) const {
//  if (!Model::CONFIG->seasonal_info().enable) {
//    return 1;
//  }
//  const auto day_of_year = TimeHelpers::day_of_year(today);
//  const auto is_rainy_period =
//      Model::CONFIG->seasonal_info().phi[location] < Constants::DAYS_IN_YEAR() / 2.0
//          ? day_of_year >= Model::CONFIG->seasonal_info().phi[location]
//                && day_of_year <= Model::CONFIG->seasonal_info().phi[location] + Constants::DAYS_IN_YEAR() / 2.0
//          : day_of_year >= Model::CONFIG->seasonal_info().phi[location]
//                || day_of_year <= Model::CONFIG->seasonal_info().phi[location] - Constants::DAYS_IN_YEAR() / 2.0;
//
//  return (is_rainy_period)
//             ? (Model::CONFIG->seasonal_info().A[location] - Model::CONFIG->seasonal_info().min_value[location])
//                       * sin(Model::CONFIG->seasonal_info().B[location] * day_of_year
//                             + Model::CONFIG->seasonal_info().C[location])
//                   + Model::CONFIG->seasonal_info().min_value[location]
//             : Model::CONFIG->seasonal_info().min_value[location];
//}
