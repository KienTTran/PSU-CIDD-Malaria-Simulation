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
    int resistant_drug_pair_id = 0;
    int resistant_type_id = 0;
    std::vector drugs = m.resistant_drug_list[resistant_drug_pair_id].second;

    std::string genotype_m = "||||YY1||KTHFI,x||||||FNCMYRIPRPCA|2";
    std::string genotype_f = "||||NY1||TTHFI,x||||||FNCMYRIPRPYA|1";
    std::string genotype_c = "||||NY1||KTHFI,x||||||FNCMYRIPRPYA|2";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;

    std::vector<Genotype*> parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    Genotype* recombined_genotype = genotype_db.get_genotype(genotype_c,&c);

    std::tuple<bool,int,int,std::string> result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 0, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[0],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id,0,"2-2"));
    EXPECT_NE(result,std::make_tuple(true,1,0,"2-2"));
    EXPECT_NE(result,std::make_tuple(true,1,1,"2-3"));
    EXPECT_NE(result,std::make_tuple(true,1,2,"2-3"));
    EXPECT_NE(result,std::make_tuple(true,1,3,""));

    //ASAQ
    resistant_drug_pair_id = 1;
    drugs = m.resistant_drug_list[resistant_drug_pair_id].second;

    //ASAQ:2-2
    genotype_m = "||||YY1||KTHFI,x||||||FNCMYRIPRPCA|1";
    genotype_f = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||NY1||KTHFI,x||||||FNCMYRIPRPYA|1";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 0, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[0],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true, resistant_drug_pair_id, 0, "2-2"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 1, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[1],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false, resistant_drug_pair_id, 1, "2-2"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 2, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[2],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false, resistant_drug_pair_id, 2, "2-2"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 3, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[3],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true, resistant_drug_pair_id, 3, ""));

    //ASAQ:2-3
    genotype_m = "||||YY1||KTHFI,x||||||FNCMYRIPRPCA|1";
    genotype_f = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||YY1||KTHFI,x||||||FNCMYRIPRPYA|1";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 0, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[0],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, 0,"2-3"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 1, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[1],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, 1,"2-3"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 2, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[2],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, 2,"2-3"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 3, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[3],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, 3,""));

    //ASAQ:2-4
    genotype_m = "||||YY1||TTHFI,x||||||FNCMYRIPRPCA|1";
    genotype_f = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||YY1||TTHFI,x||||||FNCMYRIPRPYA|1";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 0, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[0],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, 0,"2-4"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 1, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[1],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, 1,"2-4"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 2, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[2],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, 2,"2-4"));
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, 3, true);
    VLOG(0) << fmt::format("Checking: {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[3],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, 3,""));

    //Test
    //DHAPPQ
    VLOG(0) << "\n\nTest mos DHAPPQ:2-2, count DHAPPQ only";
    genotype_m = "||||NY1||TTHFI,x||||||FNCMYRIPRPCA|2";
    genotype_f = "||||NY1||TTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||NY1||TTHFI,x||||||FNCMYRIPRPYA|2";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    resistant_drug_pair_id = 0;
    resistant_type_id = 0;
    drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
    VLOG(0) << fmt::format("Checking {}: {} {} {} {}",
                           m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, resistant_type_id,"2-2"));
    for (int resistant_drug_pair_id = 0; resistant_drug_pair_id < m.resistant_drug_list.size(); resistant_drug_pair_id++) {
        auto drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
        auto resistant_types = m.resistant_drug_list[resistant_drug_pair_id].first.size();
        for (int resistant_type_id = 0; resistant_type_id < resistant_types; resistant_type_id++) {
            result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
            VLOG(0) << fmt::format("Checking {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                   std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
        }
    }

    //DHAPPQ
    VLOG(0) << "\n\nTest mos ASAQ:2-2, same genotype as above, not count DHAPPQ when checking ASAQ";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    resistant_drug_pair_id = 1;
    resistant_type_id = 0;
    drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
    VLOG(0) << fmt::format("Checking {}: {} {} {} {}",
                           m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, resistant_type_id,"0-0"));

    //DHAPPQ
    VLOG(0) << "\n\nTest mos DHAPPQ:2-2 (2), new genotype, not count DHAPPQ 1, count ASAQ 1";
    genotype_m = "||||NY1||KTHFI,x||||||FNCMYRIPRPCA|1";
    genotype_f = "||||NF1||KTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||NY1||KTHFI,x||||||FNCMYRIPRPYA|1";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    resistant_drug_pair_id = 0;
    resistant_type_id = 0;
    drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
    VLOG(0) << fmt::format("Checking {}: {} {} {} {}",
                           m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(false,resistant_drug_pair_id, resistant_type_id,"0-0"));
    resistant_drug_pair_id = 1;
    resistant_type_id = 0;
    drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
    result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
    VLOG(0) << fmt::format("Checking {}: {} {} {} {}",
                           m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                           std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
    EXPECT_EQ(result,std::make_tuple(true,resistant_drug_pair_id, resistant_type_id,"2-2"));
    for (int resistant_drug_pair_id = 0; resistant_drug_pair_id < m.resistant_drug_list.size(); resistant_drug_pair_id++) {
        auto drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
        auto resistant_types = m.resistant_drug_list[resistant_drug_pair_id].first.size();
        for (int resistant_type_id = 0; resistant_type_id < resistant_types; resistant_type_id++) {
            result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
            VLOG(0) << fmt::format("Checking {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                   std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
        }
    }

    //Test arbitrary genotype
    VLOG(0) << "\n\nTest mos two condition count arbitrary";
    genotype_m = "||||NY1||TTHFI,x||||||FNCMYRIPRPCA|2";
    genotype_f = "||||YY1||TTHFI,x||||||FNCMYRIPRPYA|1";
    genotype_c = "||||NY1||TTHFI,x||||||FNCMYRIPRPYA|2";
    VLOG(0) << genotype_m;
    VLOG(0) << genotype_f;
    VLOG(0) << genotype_c;
    parent_genotypes = {genotype_db.get_genotype(genotype_m,&c), genotype_db.get_genotype(genotype_f,&c)};
    recombined_genotype = genotype_db.get_genotype(genotype_c,&c);
    for (int resistant_drug_pair_id = 0; resistant_drug_pair_id < m.resistant_drug_list.size(); resistant_drug_pair_id++) {
        auto drugs = m.resistant_drug_list[resistant_drug_pair_id].second;
        auto resistant_types = m.resistant_drug_list[resistant_drug_pair_id].first.size();
        for (int resistant_type_id = 0; resistant_type_id < resistant_types; resistant_type_id++) {
            result = m.count_resistant_genotypes(&c, 0, parent_genotypes, recombined_genotype, drugs, resistant_drug_pair_id, resistant_type_id, true);
            VLOG(0) << fmt::format("Checking {}: {} {} {} {}",m.resistant_drug_list[resistant_drug_pair_id].first[resistant_type_id],
                                   std::get<0>(result),std::get<1>(result),std::get<2>(result),std::get<3>(result));
        }
    }

}
