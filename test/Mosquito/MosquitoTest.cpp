//
// Created by nguyentd on 3/17/2022.
//

#include <gtest/gtest.h>

#include "../MockRandom.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Mosquito/Mosquito.h"
#include "Population/Person.h"

using ::testing::_;
using ::testing::Return;

class MosquitoTest : public ::testing::Test {
protected:
  void SetUp() override {}

  // void TearDown() override {}
};

TEST_F(MosquitoTest, Initialize) {
  Mosquito m;
  Config c;
  c.read_from_file("input.yml");

  m.initialize(&c);

  EXPECT_EQ(m.genotypes_table.size(), c.number_of_tracking_days());
  EXPECT_EQ(m.genotypes_table[0].size(), c.number_of_locations());
  EXPECT_EQ(m.genotypes_table[0][0].size(), c.mosquito_config().prmc_size);
  for (Genotype* g : m.genotypes_table[0][0]) {
    EXPECT_EQ(g, nullptr);
  }

  EXPECT_EQ(c.mosquito_config().interrupted_feeding_rate.size(), c.number_of_locations());

  for (double ifr : c.mosquito_config().interrupted_feeding_rate) {
    EXPECT_EQ(ifr, 0.19);
  }
}

TEST_F(MosquitoTest, InitializeWithConfigHasCorrectNumberIFR) {
  Mosquito m;
  Config c;
  c.read_from_file("input.yml");
  c.mosquito_config().interrupted_feeding_rate = std::vector<double>(c.number_of_locations(), 0.20);
  m.initialize(&c);

  EXPECT_EQ(m.genotypes_table.size(), c.number_of_tracking_days());
  EXPECT_EQ(m.genotypes_table[0].size(), c.number_of_locations());
  EXPECT_EQ(m.genotypes_table[0][0].size(), c.mosquito_config().prmc_size);
  for (Genotype* g : m.genotypes_table[0][0]) {
    EXPECT_EQ(g, nullptr);
  }

  EXPECT_EQ(c.mosquito_config().interrupted_feeding_rate.size(), c.number_of_locations());

  for (double ifr : c.mosquito_config().interrupted_feeding_rate) {
    EXPECT_EQ(ifr, 0.20);
  }
}

TEST_F(MosquitoTest, InitializeWithConfigHasIncorrectNumberIFR) {
  Mosquito m;
  Config c;
  c.read_from_file("input.yml");
  c.mosquito_config().interrupted_feeding_rate = std::vector<double>(c.number_of_locations() + 1, 0.20);

  ASSERT_DEATH({ m.initialize(&c); }, ".*");
}

TEST_F(MosquitoTest, BuildInterruptedFeedingIndices) {
  Mosquito m;
  Config c;
  c.read_from_file("input.yml");

  //  Random r;
  //  r.initialize(1);

  MockRandom random;
  random.initialize(3);
  int number_of_expected_interrupted_feeding { 10 };
  EXPECT_CALL(random, random_poisson(_)).WillOnce(Return(number_of_expected_interrupted_feeding));

  m.initialize(&c);

  auto if_indices = m.build_interrupted_feeding_indices(&random, c.mosquito_config().interrupted_feeding_rate[0],
                                                        c.mosquito_config().prmc_size);

  ASSERT_EQ(if_indices.size(), c.mosquito_config().prmc_size);

  int sum = 0;
  for (const auto& e : if_indices) {
    sum += e;
  }

  ASSERT_EQ(sum, number_of_expected_interrupted_feeding);

  for (const auto& e : if_indices) {
    std::cout << fmt::format("{}--", e);
  }
  std::cout << std::endl;
}

TEST_F(MosquitoTest, PrmcSample) {
  Mosquito m;
  Random r;
  r.initialize(1);

  int prmc_size { 10 };
  int pop_size { 100 };

  // repeat the process 10 times
  for (int n = 0; n < 10; ++n) {
    std::vector<double> foi_distribution(pop_size, 1.0);

    // even id person will have no selection
    for (int i = 0; i < pop_size; i += 2) {
      foi_distribution[i] = 0;
    }

    std::vector<Person*> all_person;
    for (int i = 0; i < pop_size; ++i) {
      auto* p = new Person();
      p->set_last_therapy_id(i);
      all_person.push_back(p);
    }

    auto samples = r.multinomial_sampling<Person>(prmc_size, foi_distribution, all_person, true);

    EXPECT_EQ(samples.size(), prmc_size);

    for (int i = 0; i < prmc_size; ++i) {
      EXPECT_EQ(samples[i]->last_therapy_id() % 2, 1)
          << fmt::format("failed with p_id: {}", samples[i]->last_therapy_id());
    }

    for (auto* p : all_person) {
      delete p;
    }

    all_person.clear();
  }
}

TEST_F(MosquitoTest, PrmcSampleWhenThereIsNoFOI) {
  Mosquito m;
  Random r;
  r.initialize(1);

  int prmc_size { 10 };
  int pop_size { 100 };

  // repeat the process 10 times
  for (int n = 0; n < 1; ++n) {
    std::vector<double> foi_distribution(pop_size, 0);

    std::vector<Person*> all_person;
    for (int i = 0; i < pop_size; ++i) {
      auto* p = new Person();
      p->set_last_therapy_id(i);
      all_person.push_back(p);
    }

    auto samples = r.multinomial_sampling<Person>(prmc_size, foi_distribution, all_person, true);
    EXPECT_EQ(samples.size(), prmc_size);

    for (int i = 0; i < prmc_size; ++i) {
      EXPECT_EQ(samples[i], nullptr);
    }

    for (auto* p : all_person) {
      delete p;
    }

    all_person.clear();
  }
}

TEST_F(MosquitoTest, CountResistantGenotypes) {
    Mosquito m;
    Config c;
    GenotypeDatabase genotype_db;
    c.read_from_file("input.yml");

    //DHAPPQ:2-2
    int therapy_id = 8;
    auto *sc_therapy = dynamic_cast<SCTherapy *>(c.therapy_db()[therapy_id]);

    std::string g1_aas = "||||YY1||KTHFI,x||||||FNCMYRIPRPCA|2";
    std::string g2_aas = "||||NY1||TTHFI,x||||||FNCMYRIPRPYA|1";
    std::string g3_aas = "||||NY1||KTHFI,x||||||FNCMYRIPRPYA|2";

    std::vector<Genotype*> parent_genotypes = {genotype_db.get_genotype(g1_aas,&c), genotype_db.get_genotype(g2_aas,&c)};
    Genotype* g3 = genotype_db.get_genotype(g3_aas,&c);

    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "DHA-PPQ:2-2", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-4", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-4", therapy_id, true),false);

    //ASAQ
    therapy_id = 7;
    sc_therapy = dynamic_cast<SCTherapy *>(c.therapy_db()[therapy_id]);

    //ASAQ:2-2
    g1_aas = "||||YF1||KTHFI,x||||||FNCMYRIPRPCA|1";
    g2_aas = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    g3_aas = "||||YF1||KTHFI,x||||||FNCMYRIPRPYA|1";

    parent_genotypes = {genotype_db.get_genotype(g1_aas,&c), genotype_db.get_genotype(g2_aas,&c)};
    g3 = genotype_db.get_genotype(g3_aas,&c);

    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "DHA-PPQ:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-4", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-4", therapy_id, true),false);

    //ASAQ:2-3
    g1_aas = "||||YY1||KTHFI,x||||||FNCMYRIPRPCA|1";
    g2_aas = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    g3_aas = "||||YY1||KTHFI,x||||||FNCMYRIPRPYA|1";

    parent_genotypes = {genotype_db.get_genotype(g1_aas,&c), genotype_db.get_genotype(g2_aas,&c)};
    g3 = genotype_db.get_genotype(g3_aas,&c);

    m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true);
    m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true);

    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "DHA-PPQ:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-4", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-4", therapy_id, true),false);

    //ASAQ:2-4
    g1_aas = "||||YY1||TTHFI,x||||||FNCMYRIPRPCA|1";
    g2_aas = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    g3_aas = "||||YY1||TTHFI,x||||||FNCMYRIPRPYA|1";

    parent_genotypes = {genotype_db.get_genotype(g1_aas,&c), genotype_db.get_genotype(g2_aas,&c)};
    g3 = genotype_db.get_genotype(g3_aas,&c);

    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "DHA-PPQ:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-4", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-4", therapy_id, true),false);

    //ASAQ:2-4, g1 <-> g2
    g2_aas = "||||YY1||TTHFI,x||||||FNCMYRIPRPCA|1";
    g1_aas = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    g3_aas = "||||YY1||TTHFI,x||||||FNCMYRIPRPYA|1";

    parent_genotypes = {genotype_db.get_genotype(g1_aas,&c), genotype_db.get_genotype(g2_aas,&c)};
    g3 = genotype_db.get_genotype(g3_aas,&c);

    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "DHA-PPQ:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-2", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-3", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "ASAQ:2-4", therapy_id, true),true);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-2", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-3", therapy_id, true),false);
    EXPECT_EQ(m.genotype_resistant_to(&c,parent_genotypes, sc_therapy, g3, "AL:2-4", therapy_id, true),false);
}
