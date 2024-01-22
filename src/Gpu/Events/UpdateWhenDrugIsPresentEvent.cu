/* 
 * File:   UpdateByHavingDrugEvent.cpp
 * Author: Merlin
 * 
 * Created on July 31, 2013, 12:16 PM
 */

#include "UpdateWhenDrugIsPresentEvent.cuh"

#include "Gpu/Core/Scheduler.cuh"
#include "Core/Config/Config.h"
#include "Model.h"
#include "Gpu/Population/Person.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/DrugsInBlood.cuh"


GPU::UpdateWhenDrugIsPresentEvent::UpdateWhenDrugIsPresentEvent() : clinical_caused_parasite_(nullptr) {}

GPU::UpdateWhenDrugIsPresentEvent::~UpdateWhenDrugIsPresentEvent() = default;

void GPU::UpdateWhenDrugIsPresentEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p,
                                                  GPU::ClonalParasitePopulation *clinical_caused_parasite, const int &time) {
  if (scheduler!=nullptr) {

    auto *e = new UpdateWhenDrugIsPresentEvent();
    e->dispatcher = p;
    e->set_clinical_caused_parasite(clinical_caused_parasite);
    e->time = time;

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::UpdateWhenDrugIsPresentEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  if (person->drugs_in_blood()->size() > 0) {
    if (person->all_clonal_parasite_populations()->contain(clinical_caused_parasite_) && person->host_state()==
                                                                                                 GPU::Person::CLINICAL) {
      if (clinical_caused_parasite_->last_update_log10_parasite_density() <= Model::CONFIG
          ->parasite_density_level().
          log_parasite_density_asymptomatic) {
        person->set_host_state(GPU::Person::ASYMPTOMATIC);
      }
    }
    person->schedule_update_by_drug_event(clinical_caused_parasite_);
  } else {
    //        no drug in blood, reset the status of drug Update parasite
    // the drug update parasite is inserted into blood when  there is still drug in blood and also the clinical cause parasite

    for (auto i = 0; i < person->all_clonal_parasite_populations()->size(); i++) {
      const auto blood_parasite = person->all_clonal_parasite_populations()->parasites()->at(i);
      if (blood_parasite->update_function()==Model::MODEL->gpu_having_drug_update_function()) {
        //Prevent parasite density fall back to asymptomatic level early when day for parasite cleared by mono drug less than day of clinical ended.
        //This makes parasite density fall back to asymptomatic level after clinical ended.
//        person->determine_relapse_or_not(blood_parasite);
        blood_parasite->set_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
        blood_parasite->set_gpu_update_function(Model::MODEL->gpu_immunity_clearance_update_function());
      }
    }
  }
}
