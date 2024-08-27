//
// Created by nguyentd on 1/25/2022.
//

#include "Core/Random.h"
#include "gmock/gmock.h"

class MockRandom : public Random {
public:
  MOCK_METHOD(double, random_uniform, (), (override));
  MOCK_METHOD(unsigned long, random_uniform, (unsigned long range), (override));
  MOCK_METHOD(int, random_poisson, (const double& mean), (override));
  MOCK_METHOD(double, random_flat, (const double& from, const double &to), (override));
};