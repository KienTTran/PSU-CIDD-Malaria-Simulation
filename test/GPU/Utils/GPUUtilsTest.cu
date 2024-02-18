
#include <fmt/format.h>
#include "Core/Config/Config.h"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"
#include "Gpu/Utils/Utils.cuh"
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include "gtest/gtest.h"

class GPUUtilsTest : public ::testing::Test {
protected:
    std::vector<GPU::PersonUpdateInfo> h_person_update_info;
    thrust::device_vector<GPU::PersonUpdateInfo> d_person_update_info_vector;
    GPU::PersonUpdateInfo* d_person_update_info_raw_ptr = 0;
    GPU::Utils gu;
    void SetUp() override {
      gu.init();
      for(int i = 0; i < 20000; i++){
        for(int j = 0; j < 50; j++){
          GPU::PersonUpdateInfo p{};
          p.person_location = j;
          p.person_current_relative_biting_rate = 0.1*j;
          p.person_current_relative_moving_rate = 0.2*j;
          p.person_current_foi = 0.3*j;
          h_person_update_info.push_back(p);
        }
      }
      d_person_update_info_vector.resize(h_person_update_info.size());
      thrust::copy(h_person_update_info.begin(), h_person_update_info.end(), d_person_update_info_vector.begin());
      d_person_update_info_raw_ptr = thrust::raw_pointer_cast(d_person_update_info_vector.data());
    }

    void TearDown() override {
    }
};

//TEST_F(GPUUtilsTest, host_sum_bmf_by_loc_vector) {
//  printf("host_sum_bmf_by_loc_vector\n");
//  auto result = gu.host_sum_biting_moving_foi_by_loc_vector(d_person_update_info_vector);
//  for(int i = 0; i < 100; i++){
//    printf("result[%d]: loc %d sum: %f %f %f\n", i, thrust::get<0>(result[i]), thrust::get<1>(result[i]), thrust::get<2>(result[i]), thrust::get<3>(result[i]));
//  }
//}
//
//TEST_F(GPUUtilsTest, host_sum_bmf_by_loc_pointer) {
//  printf("\nhost_sum_bmf_by_loc_pointer\n");
//  d_person_update_info_raw_ptr = thrust::raw_pointer_cast(d_person_update_info_vector.data());
//  auto result = gu.host_sum_biting_moving_foi_by_loc_pointer(d_person_update_info_raw_ptr, 0, 100);
//  for(int i = 0; i < 100; i++){
//    printf("result[%d]: loc %d sum: %f %f %f\n", i, thrust::get<0>(result[i]), thrust::get<1>(result[i]), thrust::get<2>(result[i]), thrust::get<3>(result[i]));
//  }
//}

TEST_F(GPUUtilsTest, device_sum_bmf_by_loc_pointer) {
  printf("\ndevice_sum_bmf_by_loc_pointer\n");
  printf("%d\n",h_person_update_info.size());
  d_person_update_info_raw_ptr = thrust::raw_pointer_cast(d_person_update_info_vector.data());
  ThrustTVectorDevice<ThrustTuple4<int, double, double, double>> d_result = gu.device_sum_biting_moving_foi_by_loc_pointer(d_person_update_info_raw_ptr,
                                                                                                                           0,
                                                                                                                           d_person_update_info_vector.size());
  TVector<ThrustTuple4<int,double,double,double>> result(d_result.size());
  thrust::copy(d_result.begin(), d_result.end(), result.begin());
  for(int i = 0; i < 50; i++){
    printf("result[%d]: loc %d sum: %f %f %f\n", i, thrust::get<0>(result[i]), thrust::get<1>(result[i]), thrust::get<2>(result[i]), thrust::get<3>(result[i]));
  }
}