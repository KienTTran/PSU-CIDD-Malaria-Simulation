//
// Created by nguyentd on 1/27/2022.
//
#include "../MockRandom.h"
#include "Core/Config/Config.h"
#include "Parasites/Genotype.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::_;
using ::testing::Return;

TEST(GenotypePerformMutationTest, MutationNormalAAMask1) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, std::string, int, int, int, std::string>> test_cases = {
    {
        "||||000||0000000,x||||||000000000010|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        0,
        10,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPYA|1",
    },
    {
        "||||000||0000000,x||||||000000000010|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPYA|1",
        0,
        10,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||000||0000000,x||||||000000010000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        0,
        7,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRILRPCA|1",
    },
    {
        "||||000||0000000,x||||||000000010000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRILRPCA|1",
        0,
        7,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||000||1000000,x||||||000000000000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        3,
        0,
        "||||YY1||KTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||000||1000000,x||||||000000000000|0",
        "||||YY1||KTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        3,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||100||0000000,x||||||000000000000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        0,
        0,
        "||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||100||0000000,x||||||000000000000|0",
        "||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        0,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||000||0100000,x||||||000000000000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        3,
        0,
        0,
        "||||YY1||TSHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||000||0100000,x||||||0000000000100|0",
        "||||YY1||TSHFIMG,x||||||FNCMYRIPRPCA|1",
        3,
        0,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
  };

  for (const auto& [test_mask, original_str, drug_id, res_aa_id, random_aa_id, mutant_str] : test_cases) {
    c.update_mutation_mask(test_mask);
    auto origial_genotype = c.genotype_db.get_genotype(original_str, &c);

    MockRandom random;
    EXPECT_CALL(random, random_uniform(_)).WillOnce(Return(1)).WillRepeatedly(Return(1));
    EXPECT_CALL(random, random_flat(0.0,1.0)).WillOnce(Return(0.0005)).WillRepeatedly(Return(0.002));
    auto mutant_genotype = origial_genotype->perform_mutation_by_drug(&c, &random, c.drug_db()->at(drug_id),c.mutation_probability_by_locus());
    EXPECT_EQ(mutant_genotype->aa_sequence, mutant_str);
  }
}

TEST(GenotypePerformMutationTest, MutationCopyNumberVariation) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, std::string, int, int, int, std::string>> test_cases = {
    {
        "||||001||0000000,x||||||000000000000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        2,
        0,
        "||||YY2||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
    {
        "||||001||0000000,x||||||000000000000|0",
        "||||YY2||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        1,
        2,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
  };

  for (const auto& [test_mask, original_str, drug_id, res_aa_id, random_aa_id, mutant_str] : test_cases) {
    c.update_mutation_mask(test_mask);
    auto origial_genotype = c.genotype_db.get_genotype(original_str, &c);

    MockRandom random;
//    EXPECT_CALL(random, random_uniform(_)).WillOnce(Return(1)).WillRepeatedly(Return(1));
    EXPECT_CALL(random, random_flat(0.0,1.0)).WillOnce(Return(0.0005)).WillRepeatedly(Return(0.002));
    auto mutant_genotype = origial_genotype->perform_mutation_by_drug(&c, &random, c.drug_db()->at(drug_id),c.mutation_probability_by_locus());
    EXPECT_EQ(mutant_genotype->aa_sequence, mutant_str);
  }
}

TEST(GenotypePerformMutationTest, MutationNormalAAMask0) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, std::string, int, int, int, std::string>> test_cases = {
    {
        "||||000||0000000,x||||||000000000000|0",
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
        0,
        0,
        0,
        "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1",
    },
  };

  for (const auto& [test_mask, original_str, drug_id, res_aa_id, random_aa_id, mutant_str] : test_cases) {
    c.update_mutation_mask(test_mask);
    auto origial_genotype = c.genotype_db.get_genotype(original_str, &c);

    MockRandom random;
//    EXPECT_CALL(random, random_uniform(_)).WillOnce(Return(random_aa_id)).WillRepeatedly(Return(1));
//    EXPECT_CALL(random, random_flat(0.0,1.0)).WillOnce(Return(0.0005)).WillRepeatedly(Return(0.002));
    auto mutant_genotype = origial_genotype->perform_mutation_by_drug(&c, &random, c.drug_db()->at(drug_id),c.mutation_probability_by_locus());
    EXPECT_EQ(mutant_genotype->aa_sequence, mutant_str);
  }
}