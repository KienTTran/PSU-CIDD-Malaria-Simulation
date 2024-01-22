#ifndef CHANGETREATMENTCOVERAGEEVENT_CUH
#define CHANGETREATMENTCOVERAGEEVENT_CUH

#include "Gpu/Events/Event.cuh"
#include "Malaria/ITreatmentCoverageModel.h"

namespace GPU {
    class ChangeTreatmentCoverageEvent;
}
class GPU::ChangeTreatmentCoverageEvent : public GPU::Event {

 public:
  ITreatmentCoverageModel *treatment_coverage_model;

  explicit ChangeTreatmentCoverageEvent(ITreatmentCoverageModel *tcm);

  virtual ~ChangeTreatmentCoverageEvent();

  std::string name() override {
    return "ChangeTreatmentCoverageEvent";
  }

 private:
  void execute() override;
};

#endif // CHANGETREATMENTCOVERAGEEVENT_CUH
