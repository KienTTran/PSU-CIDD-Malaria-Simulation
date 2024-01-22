//
// Created by nguyentd on 3/30/2021.
//

#ifndef POMS_SRC_REPORTERS_NOVELDRUGREPORTER_CUH
#define POMS_SRC_REPORTERS_NOVELDRUGREPORTER_CUH


#include <sstream>
#include "Reporter.cuh"

namespace GPU {
    class NovelDrugReporter;
    class PersonIndexByLocationStateAgeClass;
}


class GPU::NovelDrugReporter : public GPU::Reporter {
DISALLOW_COPY_AND_ASSIGN(NovelDrugReporter)

DISALLOW_MOVE(NovelDrugReporter)

public:

  NovelDrugReporter() = default;

  ~NovelDrugReporter() override = default;

  void initialize() override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  void monthly_report() override;

private:
  void output_genotype_frequency_3(const int& number_of_genotypes, GPU::PersonIndexByLocationStateAgeClass* pi);

public:
  std::stringstream ss;
  const std::string group_sep = "-1111\t";
  const std::string sep = "\t";

  long cumulative_number_of_mutation_events_last_month = 0;

};


#endif //POMS_SRC_REPORTERS_NOVELDRUGREPORTER_CUH
