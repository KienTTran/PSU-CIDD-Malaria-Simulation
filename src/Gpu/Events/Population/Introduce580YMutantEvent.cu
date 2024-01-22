
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/Person.cuh"
#include "Introduce580YMutantEvent.cuh"

GPU::Introduce580YMutantEvent::Introduce580YMutantEvent(const int &location, const int &execute_at,
                                                   const double &fraction) : location_(location),
                                                                             fraction_(fraction) {
  time = execute_at;
}

GPU::Introduce580YMutantEvent::~Introduce580YMutantEvent() = default;

void GPU::Introduce580YMutantEvent::execute() {
  // TODO: rework on this


  auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();

  // get the approximate current frequency of 580Y in the population
  // and only fill up the different between input fraction and the current frequency

  double current_580Y_fraction = 0.0;
  double total_population_count = 0;
  for (int j = 0; j < Model::CONFIG->number_of_age_classes(); ++j) {
    for (GPU::Person* p :  pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j]) {
      total_population_count += p->all_clonal_parasite_populations()->size();
      for (GPU::ClonalParasitePopulation* pp : *p->all_clonal_parasite_populations()->parasites()) {
//        if (pp->genotype()->aa_structure()[2] == 1) {
//          current_580Y_fraction++;
//        }
      }
    }
  }

  current_580Y_fraction = total_population_count == 0 ? 0 : current_580Y_fraction / total_population_count;
  double target_fraction = fraction_ - current_580Y_fraction;
  if (target_fraction <= 0) {
    LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : Introduce 580Y Copy event with 0 cases";
    return;
  }
//  std::cout << target_fraction << std::endl;

  for (int j = 0; j < Model::CONFIG->number_of_age_classes(); ++j) {
    const auto number_infected_individual_in_ac =
      pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size() + pi->vPerson()[0][GPU::Person::CLINICAL][j].size();
    const auto number_of_importation_cases = Model::RANDOM->random_poisson(
      number_infected_individual_in_ac * target_fraction);
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

      //mutate all clonal populations
      for (auto* pp : *(p->all_clonal_parasite_populations()->parasites())) {
        auto* old_genotype = pp->genotype();
        auto* new_genotype = old_genotype->combine_mutation_to(2, 1);
        pp->set_genotype(new_genotype);
      }
    }
  }

  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : Introduce 580Y Copy event with fraction: "
            << target_fraction;
}
