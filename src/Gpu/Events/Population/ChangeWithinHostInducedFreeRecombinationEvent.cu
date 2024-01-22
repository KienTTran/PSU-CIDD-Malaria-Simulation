//
// Created by kient on 5/3/2022.
//

#include "ChangeWithinHostInducedFreeRecombinationEvent.cuh"
GPU::ChangeWithinHostInducedFreeRecombinationEvent::ChangeWithinHostInducedFreeRecombinationEvent(const bool& value,
                                                                                             const int& at_time)
    : value { value } {
  time = at_time;
}
void GPU::ChangeWithinHostInducedFreeRecombinationEvent::execute() {
  Model::CONFIG->within_host_induced_free_recombination() = value;
  LOG(INFO) << date::year_month_day{scheduler->calendar_date} << " : Change within host induced free recombination to " << value;
}
