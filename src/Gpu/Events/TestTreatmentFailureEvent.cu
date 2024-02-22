/* 
 * File:   TestTreatmentFailureEvent.cu
 * Author: Merlin
 * 
 * Created on July 31, 2013, 11:36 AM
 */

#include "TestTreatmentFailureEvent.cuh"
#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ClonalParasitePopulation.cuh"
#include "Gpu/Population/Person.cuh"


GPU::TestTreatmentFailureEvent::TestTreatmentFailureEvent() : clinical_caused_parasite_(nullptr), therapyId_(0) {}

GPU::TestTreatmentFailureEvent::~TestTreatmentFailureEvent() {
  if (executable && Model::GPU_DATA_COLLECTOR!=nullptr) {
    Model::GPU_DATA_COLLECTOR->number_of_treatments_with_therapy_ID()[therapyId_] -= 1;
  }
}

void GPU::TestTreatmentFailureEvent::schedule_event(GPU::Scheduler *scheduler, GPU::Person *p,
                                               GPU::ClonalParasitePopulation *clinical_caused_parasite,
                                               const int &time, const int &t_id) {
  if (scheduler==nullptr) {
    std::cout << "error null" << std::endl;
    assert(false);
  }
  if (scheduler!=nullptr) {
    auto *e = new TestTreatmentFailureEvent();
    e->dispatcher = p;
    e->set_clinical_caused_parasite(clinical_caused_parasite);
    e->time = time;
    e->set_therapyId(t_id);

    p->add(e);
    scheduler->schedule_individual_event(e);
  }
}

void GPU::TestTreatmentFailureEvent::execute() {
  auto *person = dynamic_cast<GPU::Person *>(dispatcher);
  LOG_IF(person->index() >= 1040 && person->index() <= 1045,INFO)
    << fmt::format("GPU::TestTreatmentFailureEvent::execute() {}",person->index());

  if (person->all_clonal_parasite_populations()->contain(clinical_caused_parasite())
      && clinical_caused_parasite_->last_update_log10_parasite_density() >
          Model::CONFIG->parasite_density_level().log_parasite_density_detectable) {

    Model::GPU_DATA_COLLECTOR->record_1_TF(person->location(), true);
    Model::GPU_DATA_COLLECTOR->record_1_treatment_failure_by_therapy(person->location(), person->age(), therapyId_);
  } else {
    Model::GPU_DATA_COLLECTOR->record_1_treatment_success_by_therapy(therapyId_);
  }
}
