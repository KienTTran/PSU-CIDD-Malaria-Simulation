//
// Created by Nguyen Tran on 3/16/2018.
//

#ifndef POMS_MFTDIFFERENTDISTRIBUTIONBYLOCATIONSTRATEGY_CUH
#define POMS_MFTDIFFERENTDISTRIBUTIONBYLOCATIONSTRATEGY_CUH

#include "IStrategy.cuh"
#include "Core/TypeDef.h"

namespace GPU{
    class MFTMultiLocationStrategy;
}

class GPU::MFTMultiLocationStrategy : public GPU::IStrategy {
 DISALLOW_COPY_AND_ASSIGN(MFTMultiLocationStrategy)

 DISALLOW_MOVE(MFTMultiLocationStrategy)

 public:
  std::vector<GPU::Therapy *> therapy_list;
  // DoubleVector2 distribution_by_location;
  DoubleVector2 distribution;
  DoubleVector2 start_distribution;
  DoubleVector2 peak_distribution;
  int starting_time{0};
  int peak_after{0};

  MFTMultiLocationStrategy();

  //    MFTStrategy(const MFTStrategy& orig);
  ~MFTMultiLocationStrategy() override;

  void add_therapy(GPU::Therapy *therapy) override;

  Therapy *get_therapy(GPU::Person *person) override;

  std::string to_string() const override;

  void update_end_of_time_step() override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;

};

#endif //POMS_MFTDIFFERENTDISTRIBUTIONBYLOCATIONSTRATEGY_CUH
