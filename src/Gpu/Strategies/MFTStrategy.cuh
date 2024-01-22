/* 
 * File:   MFTStrategy.h
 * Author: nguyentran
 *
 * Created on June 4, 2013, 11:09 AM
 */

#ifndef MFTSTRATEGY_CUH
#define MFTSTRATEGY_CUH

#include "IStrategy.cuh"
#include "Core/PropertyMacro.h"

class Random;

class Therapy;

namespace GPU{
    class MFTStrategy;
}

class GPU::MFTStrategy : public GPU::IStrategy {
 DISALLOW_COPY_AND_ASSIGN(MFTStrategy)

 DISALLOW_MOVE(MFTStrategy)

 public:
  std::vector<Therapy *> therapy_list;
  std::vector<double> distribution;

  MFTStrategy();

  //    MFTStrategy(const MFTStrategy& orig);
  virtual ~MFTStrategy();

  void add_therapy(Therapy *therapy) override;

  Therapy *get_therapy(GPU::Person *person) override;

  void update_end_of_time_step() override;

  std::string to_string() const override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;
};

#endif /* MFTSTRATEGY_CUH */
