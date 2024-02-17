//
// Created by nguyentd on 3/29/2022.
//
#include <fmt/format.h>

#include <chrono>
#include <numeric>

#include "Core/TypeDef.h"
#include "Core/Random.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Person.cuh"
#include "Gpu/Core/Random.cuh"
#include "gtest/gtest.h"

class GPURouletteTest : public ::testing::Test {
protected:
    ::Random r;
    GPU::Random gr;
    int n_location{10};
    int n_person{1000};
    int n_sample{200};
    int n_repeat{0};

    std::vector<GPU::Person *> all_person;
    std::vector<double> distribution;
    std::vector<double> sum_distribution;
    std::vector<double> sum_distribution_neg;

    void SetUp() override {
        Config c;
        c.gpu_config().n_threads = 1024;
        r.initialize(0);
        gr.init(n_location*n_person, 0);

        distribution.resize(n_location*n_person);
        sum_distribution.resize(n_location);
        sum_distribution_neg.resize(n_location);
        for(int loc = 0; loc < n_location; loc++){
            for (int i = 0; i < n_person; ++i) {
                auto *p = new GPU::Person();
                p->set_last_therapy_id(i);
                all_person.push_back(p);
                distribution[loc*n_person+i] = r.random_uniform();
            }
            int from_index = loc*n_person;
            int to_index = from_index + n_person;
            sum_distribution[loc] = std::accumulate(distribution.begin()+from_index, distribution.begin()+to_index, 0.0);
            sum_distribution_neg[loc] = -1.0;
        }
    }

    void TearDown() override {
        distribution.clear();
        std::vector<double>().swap(distribution);
        sum_distribution.clear();
        std::vector<double>().swap(sum_distribution);
        sum_distribution_neg.clear();
        std::vector<double>().swap(sum_distribution_neg);
        for (auto *p: all_person) {
            delete p;
            p = nullptr;
        }
    }
};

TEST_F(GPURouletteTest, Sampling_with_sum_distribution_0) {
    distribution = std::vector<double>(n_location*n_person, 0.0);
    for(int loc = 0; loc < n_location; loc++){
        int from_index = loc*n_person;
        int to_index = from_index + n_person;
        sum_distribution[loc] = std::accumulate(distribution.begin()+from_index, distribution.begin()+to_index, 0.0);
    }
    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution,false);

    EXPECT_EQ(results.size(), n_location*n_sample);
    EXPECT_EQ(results[0], nullptr);
}

TEST_F(GPURouletteTest, Sampling_with_no_sum_distribution) {
    distribution = std::vector<double>(n_location*n_person, 1.0);
    // even id person will have no selection
    sum_distribution = std::vector<double>(n_location, 0.0);
    for(int loc = 0; loc < n_location; loc++){
        for (int i = 0; i < n_person; i += 2) {
            distribution[loc*n_person+i] = 0;
        }
        sum_distribution[loc] = -1.0;
    }
    n_repeat = 10;

    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    ThrustTVectorDevice<double> d_sum_distribution_neg = sum_distribution_neg;
    for (int n = 0; n < n_repeat; ++n) {
        auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution_neg,false);

        EXPECT_EQ(results.size(), n_location*n_sample);

        for(int loc = 0; loc < n_location; loc++) {
            for (int i = 0; i < n_sample; ++i) {
//                std::cout << results[loc*n_sample+i]->last_therapy_id() << "\t";
                EXPECT_EQ(results[loc*n_sample+i]->last_therapy_id() % 2, 1)
                                    << fmt::format("failed with p_id: {}", results[loc*n_sample+i]->last_therapy_id());
            }
        }
        std::cout << std::endl;
    }
}

TEST_F(GPURouletteTest, Sampling_with_one_in_all) {
    distribution = std::vector<double>(n_location*n_person, 0.0);
    for(int loc = 0; loc < n_location; loc++) {
        distribution[loc*n_person + n_person - 1] = 1;
        int from_index = loc*n_person;
        int to_index = from_index + n_person;
        sum_distribution[loc] = std::accumulate(distribution.begin()+from_index, distribution.begin()+to_index, 0.0);
    }

    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    for (int n = 0; n < 10; ++n) {
        auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution, true);

        EXPECT_EQ(results.size(), n_location*n_sample);

        for(int loc = 0; loc < n_location; loc++) {
            for (int i = 0; i < n_sample; ++i) {
//                std::cout << results[loc*n_sample+i]->last_therapy_id() << "\t";
                EXPECT_EQ(results[loc*n_sample+i]->last_therapy_id(), n_person - 1)
                                    << fmt::format("failed with p_id: {}", results[loc*n_sample+i]->last_therapy_id());
            }
        }
        std::cout << std::endl;
    }
}

TEST_F(GPURouletteTest, Sampling_with_2_in_all) {
    distribution = std::vector<double> (n_location *n_person,0.0);
    for (int loc = 0; loc < n_location; loc++) {
        distribution[loc*n_person + n_person - 1] = 0.2;
        distribution[loc*n_person + n_person] = 1.8;
        sum_distribution[loc] = distribution[loc * n_person - 1] + distribution[loc * n_person];
    }
    std::vector<int> counts(n_location, 0);
    n_repeat = 10;

    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    for (int n = 0; n < n_repeat; ++n) {
        auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person,
                                                       d_sum_distribution, false);

        EXPECT_EQ(results.size(), n_location*n_sample);

        for (int loc = 0; loc < n_location; loc++) {
            for (int i = 0; i < n_sample; ++i) {
//                std::cout << results[loc * n_sample + i]->last_therapy_id() << "\t";
                EXPECT_TRUE(results[loc * n_sample + i]->last_therapy_id() == (0) ||
                            results[loc * n_sample + i]->last_therapy_id() == (n_person - 1))
                                    << fmt::format("failed with p_id: {}",
                                                   results[loc * n_sample + i]->last_therapy_id());
                if (results[loc * n_sample + i]->last_therapy_id() == (n_person - 1)) {
                    counts[loc]++;
                }
            }
            std::cout << std::endl;
        }
    }

    for (int loc = 0; loc < n_location; loc++) {
        std::cout << fmt::format("Expected - Actual freq: {} - {}",
                                 distribution[loc * n_person + n_person - 1] / sum_distribution[loc],
                                 counts[loc] / (double) (n_repeat * n_sample)) << std::endl;
//        EXPECT_NEAR(distribution[loc * n_person + n_person - 1] / sum_distribution[loc],
//                    counts[loc] / (double) (n_repeat * n_sample), 0.02);
    }
}

TEST_F(GPURouletteTest, Sampling_with_4_in_all) {
    distribution = std::vector<double>(n_location*n_person, 0.0);
    for(int loc = 0; loc < n_location; loc++) {
        distribution[loc*n_person] = 0.2;
        distribution[loc*n_person + 100] = 0.5;
        distribution[loc*n_person + 200] = 0.8;
        distribution[loc*n_person + 300] = 1.5;
        sum_distribution[loc] = distribution[loc*n_person] + distribution[loc*n_person + 100]
                                + distribution[loc*n_person + 200] + distribution[loc*n_person + 300];
    }
    n_repeat = 1000;
    std::map<int, int> count = {
            {0,   0},
            {100, 0},
            {200, 0},
            {300, 0},
    };
    std::vector<std::map<int, int>> counts = std::vector<std::map<int, int>>(n_location);
    for (int i = 0; i < n_location; ++i) {
        counts[i] = count;
    }

    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    for (int n = 0; n < n_repeat; ++n) {
        auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution, true);

        EXPECT_EQ(results.size(), n_location*n_sample);

        for (int loc = 0; loc < n_location; loc++) {
            for (int i = 0; i < n_sample; ++i) {
                EXPECT_TRUE(results[loc * n_sample + i]->last_therapy_id() == 0 ||
                            results[loc * n_sample + i]->last_therapy_id() == 100
                            || results[loc * n_sample + i]->last_therapy_id() == 200 ||
                            results[loc * n_sample + i]->last_therapy_id() == 300)
                                    << fmt::format("failed with p_id: {}", results[i]->last_therapy_id());
                counts[loc][results[loc * n_sample + i]->last_therapy_id()]++;
            }
        }
    }

    for(int loc = 0; loc < n_location; loc++){
        std::cout << fmt::format("{} trials, {} samples per trial", n_repeat, n_sample) << std::endl;
        for (const auto &[key, value]: counts[loc]) {
            std::cout << fmt::format("Key: {} \tExpected - Actual freq: {} - {}", key,
                                     distribution[key] / sum_distribution[loc],
                                     value / (double) (n_repeat * n_sample))
                      << std::endl;
            EXPECT_NEAR(distribution[key] / sum_distribution[loc], value / (double) (n_repeat * n_sample), 0.01);
        }
    }
}

TEST_F(GPURouletteTest, compare_with_multi_normial) {
    distribution = std::vector<double>(n_location*n_person, 0.0);
    for(int loc = 0; loc < n_location; loc++) {
        for (int i = 0; i < n_person; ++i) {
            distribution[loc*n_person+i] = r.random_uniform();
        }
        int from_index = loc*n_person;
        int to_index = from_index + n_person;
        sum_distribution[loc] = std::accumulate(distribution.begin()+from_index, distribution.begin()+to_index, 0.0);
    }
    n_repeat = 1000;

    // ==================== roulette sampling ===========================
    auto start = std::chrono::high_resolution_clock::now();

    ThrustTVectorDevice<double> d_distribution = distribution;
    ThrustTVectorDevice<double> d_sum_distribution = sum_distribution;
    for (int n = 0; n < n_repeat; ++n) {
        auto results = gr.roulette_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution, true);

        EXPECT_EQ(results.size(), n_location*n_sample);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << fmt::format("{} locations, {} trials, {} samples per trial", n_location, n_repeat, n_sample) << std::endl;
    std::cout << fmt::format("Roulette Sampling times: {}ms", duration.count()) << std::endl;

    // ==================== multinomial sampling ===========================
    start = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < n_repeat; ++n) {
        auto results = gr.multinomial_sampling<GPU::Person>(n_location, n_sample, d_distribution, all_person, d_sum_distribution, true);

        EXPECT_EQ(results.size(), n_location*n_sample);
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << fmt::format("{} locations, {} trials, {} samples per trial", n_location, n_repeat, n_sample) << std::endl;
    std::cout << fmt::format("Multinomial Sampling times: {}ms", duration.count()) << std::endl;
}

