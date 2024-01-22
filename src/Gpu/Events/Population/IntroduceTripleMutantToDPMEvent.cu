//
// Created by nguyentd on 11/3/2020.
//

#include "IntroduceTripleMutantToDPMEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Population/SingleHostClonalParasitePopulations.h"

GPU::IntroduceTrippleMutantToDPMEvent::IntroduceTrippleMutantToDPMEvent(
    const int& location, const int& execute_at,
    const double& fraction
)
    : location_(location),
      fraction_(fraction) {
  time = execute_at;
}

GPU::IntroduceTrippleMutantToDPMEvent::~IntroduceTrippleMutantToDPMEvent() = default;

void GPU::IntroduceTrippleMutantToDPMEvent::execute() {
  auto* pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>();


  for (int j = 0; j < Model::CONFIG->number_of_age_classes(); ++j) {
    const auto number_infected_individual_in_ac =
        pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size() + pi->vPerson()[0][GPU::Person::CLINICAL][j].size();
    const auto number_of_importation_cases = Model::RANDOM->random_poisson(
        number_infected_individual_in_ac * fraction_
    );
    if (number_of_importation_cases == 0) {
      continue;
    }
    for (auto i = 0; i < number_of_importation_cases; i++) {

      const size_t index = Model::RANDOM->random_uniform(number_infected_individual_in_ac);

      GPU::Person* p = nullptr;
      if (index < pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size()) {
        p = pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j][index];
      } else {
        p = pi->vPerson()[0][GPU::Person::CLINICAL][j][index - pi->vPerson()[0][GPU::Person::ASYMPTOMATIC][j].size()];
      }

      // get random parasite population
      //mutate all
      for (auto* pp : *(p->all_clonal_parasite_populations()->parasites())) {
        // TODO: rework on this
        auto* old_genotype = pp->genotype();
        auto* new_genotype = old_genotype->combine_mutation_to(2, 1)
                                         ->combine_mutation_to(3, 1);

        // mutate to double copy of PfMdr 2 copies
//        auto mdr_gene_allele_value = new_genotype->aa_structure[1];
//        if (mdr_gene_allele_value < 4) {
//          new_genotype = new_genotype->combine_mutation_to(1, mdr_gene_allele_value + 4);
//        }

        pp->set_genotype(new_genotype);
      }
    }
  }

  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : IntroduceTrippleMutantToDPMEvent";
}
