/*
 * File:   Person.h
 * Author: nguyentran
 *
 * Created on March 22, 2013, 2:25 PM
 */

#ifndef PERSON_CUH
#define PERSON_CUH

#include "Gpu/Core/Dispatcher.cuh"
#include "Core/ObjectPool.h"
#include "Core/PropertyMacro.h"
#include "Gpu/Population/Properties/PersonIndexAllHandler.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationMovingLevelHandler.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClassHandler.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Parasites/Genotype.cuh"
#include "Gpu/Population/Properties/PersonIndexGPUHandler.hxx"

namespace GPU{
    class SingleHostClonalParasitePopulations;
    class ClonalParasitePopulation;
    class ParasiteDensityUpdateFunction;
    class ImmuneSystem;
    class Genotype;
    class Person;
    class DrugsInBlood;
    class Population;
    class Scheduler;
    class Event;
    class ModelDataCollector;
    class Therapy;
    class DrugType;
    class PersonIndexGPU;
}

class Model;

class Config;

class Random;

class GPU::Person : public GPU::PersonIndexAllHandler,
               public GPU::PersonIndexByLocationStateAgeClassHandler,
               public GPU::PersonIndexByLocationMovingLevelHandler,
               public GPU::PersonIndexGPUHandler,
               public GPU::Dispatcher {
public:

  enum Property {
    LOCATION = 0,
    HOST_STATE,
    AGE,
    AGE_CLASS,
    BITTING_LEVEL,
    MOVING_LEVEL,
    EXTERNAL_POPULATION_MOVING_LEVEL
  };

  enum HostStates { SUSCEPTIBLE = 0, EXPOSED = 1, ASYMPTOMATIC = 2, CLINICAL = 3, DEAD = 4, NUMBER_OF_STATE = 5 };

  POINTER_PROPERTY(GPU::Population, population)

  PROPERTY_HEADER(int, location)

  PROPERTY_REF(int, residence_location)

  PROPERTY_HEADER(HostStates, host_state)

  PROPERTY_HEADER(int, age)

  PROPERTY_HEADER(int, age_class)

  PROPERTY_HEADER(long, id)

  PROPERTY_HEADER(int, index)

  // birthday has the unit of time in the scheduler
  // if birthday is -100 which is that person was born 100 day before the simulation start
  PROPERTY_REF(int, birthday)

  POINTER_PROPERTY_HEADER(GPU::ImmuneSystem, immune_system)

  POINTER_PROPERTY(GPU::SingleHostClonalParasitePopulations, all_clonal_parasite_populations)

  VIRTUAL_PROPERTY_REF(int, latest_update_time)

  PROPERTY_HEADER(int, moving_level)

  POINTER_PROPERTY(GPU::DrugsInBlood, drugs_in_blood)

  POINTER_PROPERTY(GPU::Genotype, liver_parasite_type)

  POINTER_PROPERTY(IntVector, today_infections)

  POINTER_PROPERTY(IntVector, today_target_locations)

  PROPERTY_REF(int, number_of_times_bitten)

  PROPERTY_REF(int, number_of_trips_taken)
  //    PROPERTY_REF(bool, is_tracking_treatment_number);
  PROPERTY_REF(int, last_therapy_id)

  PROPERTY_REF(std::vector<double>, prob_present_at_mda_by_age)

  PROPERTY_REF(int, location_col)
  PROPERTY_REF(int, location_row)
  PROPERTY_REF(glm::mat4, model)
  PROPERTY_REF(glm::vec4, color)

  GPU::PersonIndexGPU* person_index_gpu;

public:
  std::map<int, double> starting_drug_values_for_MAC;
  double innate_relative_biting_rate { 0 };
  double current_relative_biting_rate { 0 };

public:
  Person();

  //    Person(const Person& orig);
  virtual ~Person();

  void init() override;

  //    Model* model();
  //    Scheduler* scheduler();
  //    Config* config();
  //    Random* random();

  void NotifyChange(const Property &property, const void *oldValue, const void *newValue);

  virtual void increase_age_by_1_year();

  //    BloodParasite* add_new_parasite_to_blood(Genotype* parasite_type);
  GPU::ClonalParasitePopulation *add_new_parasite_to_blood(GPU::Genotype *parasite_type);

  static double relative_infectivity(const double &log10_parasite_density);

  virtual double get_probability_progress_to_clinical();

  virtual bool will_progress_to_death_when_receive_no_treatment();

  virtual bool will_progress_to_death_when_recieve_treatment();

  void cancel_all_other_progress_to_clinical_events_except(GPU::Event *event) const;

  void cancel_all_events_except(GPU::Event *event) const;
  //    void record_treatment_failure_for_test_treatment_failure_events();

  void change_all_parasite_update_function(GPU::ParasiteDensityUpdateFunction *from,
                                           GPU::ParasiteDensityUpdateFunction *to);

  int complied_dosing_days(const int &dosing_day) const;

  void receive_therapy(GPU::Therapy *therapy, GPU::ClonalParasitePopulation *clinical_caused_parasite,
                       bool is_part_of_MAC_therapy = false);

  void add_drug_to_blood(GPU::DrugType *dt, const int &dosing_days, bool is_part_of_MAC_therapy = false);

  void schedule_progress_to_clinical_event_by(GPU::ClonalParasitePopulation *blood_parasite);

  void schedule_test_treatment_failure_event(GPU::ClonalParasitePopulation *blood_parasite, const int &testing_day,
                                             const int &t_id = 0);

  void schedule_update_by_drug_event(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void schedule_end_clinical_event(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void schedule_end_clinical_by_no_treatment_event(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void schedule_relapse_event(GPU::ClonalParasitePopulation *clinical_caused_parasite, const int &time_until_relapse);

  void schedule_move_parasite_to_blood(GPU::Genotype *genotype, const int &time);

  void schedule_mature_gametocyte_event(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void change_state_when_no_parasite_in_blood();

  void determine_relapse_or_not(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void determine_clinical_or_not(GPU::ClonalParasitePopulation *clinical_caused_parasite);

  void update() override;

  void update_current_state();

  void randomly_choose_parasite();

  void infected_by(const int &parasite_type_id);

  void randomly_choose_target_location();

  void schedule_move_to_target_location_next_day_event(const int &location);

  bool has_return_to_residence_event() const;

  void cancel_all_return_to_residence_events() const;

  bool has_detectable_parasite() const;

  void increase_number_of_times_bitten();

  void move_to_population(GPU::Population *target_population);

  bool has_birthday_event() const;

  bool has_update_by_having_drug_event() const;

  double get_age_dependent_biting_factor() const;

  void update_relative_bitting_rate();

  double p_infection_from_an_infectious_bite() const;

  bool isGametocytaemic() const;

  void generate_prob_present_at_mda_by_age();

  double prob_present_at_mda();

  bool has_effective_drug_in_blood() const;

  static double draw_random_relative_biting_rate(::Random *pRandom, Config *pConfig);

    void generate_render_entity(int location, bool is_circulate = false);
};

#endif /* PERSON_CUH */
