//
// Created by kient on 8/22/2023.
//

#include "ChangeMutationProbabilityByLocusEvent.cuh"

GPU::ChangeMutationProbabilityByLocusEvent::ChangeMutationProbabilityByLocusEvent(const double& value,
                                                                                             const int& at_time)
        : value { value } {
    time = at_time;
}
void GPU::ChangeMutationProbabilityByLocusEvent::execute() {
    Model::CONFIG->mutation_probability_by_locus() = value;
    LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : Change mutation probability by locus to " << value;
}