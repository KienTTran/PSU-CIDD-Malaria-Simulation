
#include "Core/Config/Config.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::ElementsAre;

TEST(ConfigTest, ReadGPUConfigTest) {
    Config c;
    c.read_from_file("input.yml");

    EXPECT_EQ(c.gpu_config().n_threads, 1024);
    EXPECT_EQ(c.render_config().width, 1280);
    EXPECT_EQ(c.render_config().height, 720);
    EXPECT_EQ(c.debug_config().log_interval, 100);
}