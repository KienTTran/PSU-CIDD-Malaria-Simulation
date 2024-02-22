/* 
 * File:   EndClinicalDueToDrugResistanceEvent.cu
 * Author: Merlin
 * 
 * Created on July 31, 2013, 11:24 AM
 */

#include "EndClinicalDueToDrugResistanceEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/Person.cuh"

GPU::EndClinicalDueToDrugResistanceEvent::EndClinicalDueToDrugResistanceEvent() : clinical_caused_parasite_(nullptr) {}

GPU::EndClinicalDueToDrugResistanceEvent::~EndClinicalDueToDrugResistanceEvent() = default;

void GPU::EndClinicalDueToDrugResistanceEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p,
                                                         GPU::ClonalParasitePopulation *clinical_caused_parasite,
                                                         const int &time) {
  if (scheduler!=nullptr) {
    auto *e = new GPU::EndClinicalDueToDrugResistanceEvent();
    e->dispatcher = p;
    e->set_clinical_caused_parasite(clinical_caused_parasite);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::EndClinicalDueToDrugResistanceEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  if (person->all_clonal_parasite_populations()->size()==0) {
    person->change_state_when_no_parasite_in_blood();

  } else {
    //still have parasite in blood
    person->immune_system()->set_increase(true);
    person->set_host_state(GPU::Person::ASYMPTOMATIC);

    if (person->all_clonal_parasite_populations()->contain(clinical_caused_parasite_)) {
      clinical_caused_parasite_->set_last_update_log10_parasite_density(
          Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic);

      person->determine_relapse_or_not(clinical_caused_parasite_);

    }

    //        person->determine_relapse_or_not(clinical_caused_parasite_);
    //        person->determine_clinical_or_not(clinical_caused_parasite_);

    //        if (clinical_caused_parasite_->last_update_log10_parasite_density() < Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic) {
    //            std::cout << clinical_caused_parasite_->last_update_log10_parasite_density()<< std::endl;
    //            assert(false);
    //        }
    //        clinical_caused_parasite_->set_last_update_log10_parasite_density(Model::CONFIG->parasite_density_level().log_parasite_density_asymptomatic);
    //        clinical_caused_parasite_->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
    //        //        std::cout << "hello" << std::endl;
  }

}
