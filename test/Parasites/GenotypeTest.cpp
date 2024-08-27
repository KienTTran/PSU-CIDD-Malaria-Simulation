//
// Created by nguyentd on 12/7/2021.
//

#include <gtest/gtest.h>

#include "Core/Config/Config.h"
#include "Parasites/Genotype.h"

class GenotypeTest : public ::testing::Test {
protected:
  void SetUp() override {}

  // void TearDown() override {}
};

TEST_F(GenotypeTest, InitializeWithAASequence) {
  Genotype g("||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1");

  EXPECT_EQ(g.pf_genotype_str[0].size(), 0);
  EXPECT_EQ(g.pf_genotype_str[4][0], "NY1");
  EXPECT_EQ(g.pf_genotype_str[6][0], "KTHFIMG");
  EXPECT_EQ(g.pf_genotype_str[6][1], "x");
  EXPECT_EQ(g.pf_genotype_str[12][0], "FNCMYRIPRPCA");
  EXPECT_EQ(g.pf_genotype_str[13][0], "1");

  EXPECT_EQ(g.get_aa_sequence(), "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1");
}

TEST_F(GenotypeTest, CheckValid) {
  Config c;
  c.read_from_file("input.yml");

  Genotype g("||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  EXPECT_TRUE(g.is_valid(c.pf_genotype_info()));

  Genotype g2("||||NY2||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  EXPECT_TRUE(g2.is_valid(c.pf_genotype_info()));

  Genotype g3("||||NY2||KTHFIMG,x||||||FNCMYRIPRPCA|2");
  EXPECT_TRUE(g3.is_valid(c.pf_genotype_info()));
}

TEST_F(GenotypeTest, CheckInvalid) {
  Config c;
  c.read_from_file("input.yml");

  Genotype g1("||||NY1||KTHFI||||||FNCMYRIPRPCA|1");
  EXPECT_FALSE(g1.is_valid(c.pf_genotype_info()));

  Genotype g2("|");
  EXPECT_FALSE(g2.is_valid(c.pf_genotype_info()));

  Genotype g3("||||MY1||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  EXPECT_FALSE(g3.is_valid(c.pf_genotype_info()));

  Genotype g4("||||NY3||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  EXPECT_FALSE(g4.is_valid(c.pf_genotype_info()));

  Genotype g5("||||NYC2||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  EXPECT_FALSE(g5.is_valid(c.pf_genotype_info()));

  Genotype g6("||||NY2||KTHFIMG,x||||||FNCMYRIPRPC|Y1");
  EXPECT_FALSE(g6.is_valid(c.pf_genotype_info()));
}

TEST_F(GenotypeTest, CalculateCostOfResistance) {
  Config c;
  c.read_from_file("input.yml");

  Genotype g("||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1");
  g.calculate_daily_fitness(c.pf_genotype_info());
  EXPECT_DOUBLE_EQ(g.daily_fitness_multiple_infection, 1);

  Genotype g1("||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1");
  g1.calculate_daily_fitness(c.pf_genotype_info());
  EXPECT_NEAR(g1.daily_fitness_multiple_infection, pow(1 - 0.0005, 1),0.0001);

  Genotype g2("||||YF2||TTHFIMG,x||||||FNCMYRIPRPCA|1");
  g2.calculate_daily_fitness(c.pf_genotype_info());
  EXPECT_NEAR(g2.daily_fitness_multiple_infection, pow(1 - 0.0005, 3) * (1 - 0.0005),0.0001);
}

TEST_F(GenotypeTest, Calculate_Base_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_samples = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 0, pow(0.75, 25) },
    { "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.6, 20) },
    { "||||NF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.5, 19) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 3, pow(0.58, 15) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 4, pow(0.45, 15) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 5, pow(1.08, 15) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.72, 19) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 7, pow(0.55, 15) }
  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_samples) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.0000000001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_AS_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 0, pow(0.75, 25) },
    // single as mutation
    { "||||NY1||KTHFIMG,x||||||INCMYRIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FYCMYRIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNYMYRIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCIYRIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMHRIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYTIPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRTPRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRILRPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPHPCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRLCA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPYA|1", 0, pow(1.2, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCV|1", 0, pow(1.2, 25) },

    // double mutation
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPYV|1", 0, pow(1.32, 25) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPHLCA|1", 0, pow(1.32, 25) },

    // triple mutation
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRLYV|1", 0, pow(1.452, 25) },

    // triple mutation on k13 and none relevant genes
    { "||||YF2||KTHFIMG,x||||||FNCMYRIPRLYV|1", 0, pow(1.452, 25) },
    { "||||YF2||KTHFIMG,x||||||FNCMYRIPRLYV|2", 0, pow(1.452, 25) },

    // quadruple mutation on k13 and none relevant genes
    { "||||YF2||KTHFIMG,x||||||FNCMYRIPHLYV|2", 0, pow(1.5972, 25) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.000001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_Piperaquine_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 3, pow(0.58, 15) },
    // mutation on none relevant genes
    { "||||YF2||TTHFIMG,x||||||FNCMYRIPRPYA|1", 3, pow(0.58, 15) },

    // copy number variation
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|2", 3, pow(0.7946, 15) },

    // copy number variation and none relevant genes
    { "||||NY2||KTHFIMG,x||||||FNCMYRIPRPYA|2", 3, pow(0.7946, 15) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_Mefloquine_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 4, pow(0.45, 15) },
    // mutation on none relevant genes
    { "||||YF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 4, pow(0.45, 15) },

    // copy number variation
    { "||||NY2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 4, pow(1.1, 15) },

    // copy number variation and none relevant genes
    { "||||YF2||TTHFIMG,x||||||FNCMYRIPRPYA|1", 4, pow(1.1, 15) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_SP_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 5, pow(1.08, 15) },
    // mutations has no effect
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|2", 5, pow(1.08, 15) },
    { "||||YF2||TTHFIMG,x||||||FNCMYRIPRPYA|1", 5, pow(1.08, 15) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_CQ_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.72, 19) },

    // only N85Y in mdr1 has effect
    { "||||NF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.72, 19) },
    { "||||NY2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.72, 19) },
    { "||||NF2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.72, 19) },
    { "||||YY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.9, 19) },
    { "||||YF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.9, 19) },
    { "||||YF2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(0.9, 19) },

    // K76T
    { "||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.152, 19) },
    { "||||NY1||TSHFIMG,x||||||FNCMYRIPRPYA|2", 6, pow(1.152, 19) },
    { "||||NY2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.152, 19) },
    { "||||NF1||TSHFIMG,x||||||FNCMYRIPRPYA|2", 6, pow(1.152, 19) },
    { "||||NF2||TSHFIMG,x||||||FNCMYRIPRPYA|2", 6, pow(1.152, 19) },

    // 85Y + 76T
    { "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.44, 19) },
    { "||||YF1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.44, 19) },
    { "||||YY2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.44, 19) },
    { "||||YF2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 6, pow(1.44, 19) },

    //    1.19 # TODO: test this value
  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_Lumefantrine_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    // base ec50
    { "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.6, 20) },

    // mutation on non-relevant genes
    { "||||YY1||TTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(0.6, 20) },

    // single mutation N86
    { "||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.75, 20) },
    { "||||NY1||TTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(0.75, 20) },

    // single mutation 184F
    { "||||YF1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.75, 20) },
    { "||||YF1||TTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(0.75, 20) },

    // single mutation K76
    { "||||YY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.66, 20) },
    { "||||YY1||KTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(0.66, 20) },

    // double copy YYYY + 76T
    { "||||YY2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.78, 20) },

    // double copy YYYY + K76
    { "||||YY2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.858, 20) },

    // N86 Y184 K76
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.825, 20) },
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(0.825, 20) },

    // NYNY - K76
    { "||||NY2||KTHFIMG,x||||||FNCMYRIPRPYA|2", 1, pow(1.0725, 20) },

    // N86 - 184F 76T
    { "||||NF1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.7875, 20) },

    // NFNF 76T
    { "||||NF2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(1.02375, 20) },

    // NF-NF K76
    { "||||NF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(0.86625, 20) },
    { "||||NF2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 1, pow(1.126125, 20) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}

TEST_F(GenotypeTest, Calculate_Amodiaquine_EC50) {
  Config c;
  c.read_from_file("input.yml");

  std::vector<std::tuple<std::string, int, double>> test_cases = {
    { "||||NY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.6, 19) },
    { "||||YY1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.852, 19) },
    { "||||NF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.5, 19) },
    { "||||YF1||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.71, 19) },

    { "||||NY2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.6, 19) },
    { "||||YY2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.852, 19) },
    { "||||NF2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.5, 19) },
    { "||||YF2||KTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.71, 19) },

    { "||||NY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.72, 19) },
    { "||||YY1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(1.0224, 19) },
    { "||||NF1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.6, 19) },
    { "||||YF1||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.852, 19) },

    { "||||NY2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.72, 19) },
    { "||||YY2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(1.0224, 19) },
    { "||||NF2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.6, 19) },
    { "||||YF2||TTHFIMG,x||||||FNCMYRIPRPCA|1", 2, pow(0.852, 19) },

  };

  for (const auto& [gene_str, drug_id, ec50_p_n] : test_cases) {
    Genotype g(gene_str);
    g.calculate_EC50_power_n(c.pf_genotype_info(), c.drug_db());
    EXPECT_NEAR(g.get_EC50_power_n(c.drug_db()->at(drug_id)), ec50_p_n, 0.00001)
        << fmt::format("{}-{}-{}", gene_str, drug_id, ec50_p_n);
  }
}
