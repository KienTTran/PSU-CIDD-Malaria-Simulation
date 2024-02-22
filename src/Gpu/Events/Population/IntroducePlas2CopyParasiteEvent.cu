//
// Created by Nguyen Tran on 2019-05-01.
//

#include "IntroducePlas2CopyParasiteEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Population/SingleHostClonalParasitePopulations.h"

GPU::IntroducePlas2CopyParasiteEvent::IntroducePlas2CopyParasiteEvent(const int &location, const int &execute_at,
                                                                 const double &fraction)
  : location_(location),
    fraction_(fraction) {
  time = execute_at;
}

GPU::IntroducePlas2CopyParasiteEvent::~IntroducePlas2CopyParasiteEvent() = default;

void GPU::IntroducePlas2CopyParasiteEvent::execute() {
  auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();


  for (int j = 0; j < Model::CONFIG->number_of_age_classes(); ++j) {
    const auto number_infected_individual_in_ac =
      pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size() + pi->vPerson()[0][GPU::Person::CLINICAL][j].size();
    const auto number_of_importation_cases = Model::RANDOM->random_poisson(
      number_infected_individual_in_ac * fraction_);
    if (number_of_importation_cases == 0)
      continue;
    for (auto i = 0; i < number_of_importation_cases; i++) {

      const size_t index = Model::RANDOM->random_uniform(number_infected_individual_in_ac);

      GPU::Person* p = nullptr;
      if (index < pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size()) {
        p = pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j][index];
      } else {
        p = pi->vPerson()[0][GPU::Person::CLINICAL][j][index - pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size()];
      }

      // get random parasite population
//mutate one
//    auto pp_index = Model::RANDOM->random_uniform(p->all_clonal_parasite_populations()->size());
//    auto* pp = p->all_clonal_parasite_populations()->parasites()->at(pp_index);
//
//    auto* old_genotype = pp->genotype();
//    auto* new_genotype = old_genotype->combine_mutation_to(3, 1);
//    pp->set_genotype(new_genotype);

      //mutate all
      for (auto* pp : *(p->all_clonal_parasite_populations()->parasites())) {
        auto* old_genotype = pp->genotype();
        auto* new_genotype = old_genotype->combine_mutation_to(3, 1);
        pp->set_genotype(new_genotype);
      }
    }
  }

  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : Introduce Plas2 Copy";
}