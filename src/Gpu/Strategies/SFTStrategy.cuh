/* 
 * File:   SFTStrategy.h
 * Author: nguyentran
 *
 * Created on June 3, 2013, 8:00 PM
 */

#ifndef SFTSTRATEGY_CUH
#define SFTSTRATEGY_CUH

#include "Core/PropertyMacro.h"
#include "IStrategy.cuh"

namespace GPU {
    class SFTStrategy;
}

class GPU::SFTStrategy : public GPU::IStrategy {
 DISALLOW_COPY_AND_ASSIGN(SFTStrategy)

 DISALLOW_MOVE(SFTStrategy)

 VIRTUAL_PROPERTY_REF(std::vector<GPU::Therapy *>, therapy_list)

 public:
  SFTStrategy();

  //    SFTStrategy(const SFTStrategy& orig);
  virtual ~SFTStrategy();

  virtual std::vector<GPU::Therapy *> &get_therapy_list();

  void add_therapy(GPU::Therapy *therapy) override;

  Therapy *get_therapy(GPU::Person *person) override;

  std::string to_string() const override;

  void update_end_of_time_step() override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;

};

#endif /* SFTSTRATEGY_CUH */
