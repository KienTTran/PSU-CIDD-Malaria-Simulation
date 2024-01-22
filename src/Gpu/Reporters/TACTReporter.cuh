//
// Created by nguyentd on 7/8/2020.
//

#ifndef POMS_TACTREPORTER_CUH
#define POMS_TACTREPORTER_CUH


#include <Core/PropertyMacro.h>
#include <sstream>
#include "Reporter.cuh"

namespace GPU {
    class TACTReporter;
    class PersonIndexByLocationStateAgeClass; }

class GPU::TACTReporter : public GPU::Reporter {
DISALLOW_COPY_AND_ASSIGN(TACTReporter)

DISALLOW_MOVE(TACTReporter)

public:

  TACTReporter() = default;

  ~TACTReporter() override = default;

  void initialize() override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  void monthly_report() override;

private:
  void output_genotype_frequency_3(
      const int& number_of_genotypes,
      GPU::PersonIndexByLocationStateAgeClass* pi
  );

public:
  std::stringstream ss;
  const std::string group_sep = "-1111\t";
  const std::string sep = "\t";

  long cumulative_number_of_mutation_events_last_month = 0;

};


#endif //POMS_TACTREPORTER_CUH
