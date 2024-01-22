
#ifndef MMCREPORTER_CUH
#define MMCREPORTER_CUH

#include "Reporter.cuh"
#include <sstream>

namespace GPU {
    class MMCReporter;
    class PersonIndexByLocationStateAgeClass;
}

class GPU::MMCReporter : public GPU::Reporter {
DISALLOW_COPY_AND_ASSIGN(MMCReporter)

DISALLOW_MOVE(MMCReporter)

public:
  std::stringstream ss;
  const std::string group_sep = "-1111\t";
  const std::string sep = "\t";

  MMCReporter();

  ~MMCReporter() override = default;

  void initialize() override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  void print_treatment_failure_rate_by_therapy();

  void print_ntf_by_location();

  void monthly_report() override;

  void print_EIR_PfPR_by_location();

};

#endif // MMCREPORTER_CUH
