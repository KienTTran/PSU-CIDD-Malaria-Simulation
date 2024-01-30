#ifndef NESTEDMFTSTRATEGY_CUH
#define NESTEDMFTSTRATEGY_CUH

#include "IStrategy.cuh"

namespace GPU{
    class NestedMFTStrategy;
}

class GPU::NestedMFTStrategy : public GPU::IStrategy {
 DISALLOW_COPY_AND_ASSIGN(NestedMFTStrategy)

 DISALLOW_MOVE(NestedMFTStrategy)

 public:
  std::vector<GPU::IStrategy *> strategy_list;
  std::vector<double> distribution;
  std::vector<double> start_distribution;
  std::vector<double> peak_distribution;
  int starting_time{0};
  int peak_after{0};

  NestedMFTStrategy() : GPU::IStrategy("NestedMFTStrategy", NestedMFT) {}

  virtual ~NestedMFTStrategy() = default;

  virtual void add_strategy(GPU::IStrategy *strategy);

  void add_therapy(GPU::Therapy *therapy) override;

  Therapy *get_therapy(GPU::Person *person) override;

  std::string to_string() const override;

  void adjust_started_time_point(const int &current_time) override;

  void update_end_of_time_step() override;

  void monthly_update() override;

  void adjust_distribution(const int &time);
};

#endif // NESTEDMFTSTRATEGY_CUH
