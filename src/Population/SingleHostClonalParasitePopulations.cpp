/*
 * File:   ParasitePopulation.cpp
 * Author: Merlin
 *
 * Created on July 11, 2013, 1:53 PM
 */

#include "SingleHostClonalParasitePopulations.h"

#include <cmath>

#include "ClonalParasitePopulation.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "DrugsInBlood.h"
#include "Helpers/NumberHelpers.h"
#include "Helpers/ObjectHelpers.h"
#include "MDC/ModelDataCollector.h"
#include "Model.h"
#include "Parasites/Genotype.h"
#include "Person.h"
#include "Therapies/Drug.h"

OBJECTPOOL_IMPL(SingleHostClonalParasitePopulations)

SingleHostClonalParasitePopulations::SingleHostClonalParasitePopulations(Person* person)
    : person_(person),
      parasites_(nullptr),
      log10_total_relative_density_(ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY) {}

void SingleHostClonalParasitePopulations::init() {
  parasites_ = new std::vector<ClonalParasitePopulation*>();
}

SingleHostClonalParasitePopulations::~SingleHostClonalParasitePopulations() {
  if (parasites_ != nullptr) {
    clear();
    ObjectHelpers::delete_pointer<std::vector<ClonalParasitePopulation*>>(parasites_);
  }

  person_ = nullptr;
}

void SingleHostClonalParasitePopulations::clear() {
  if (parasites_->empty()) {
    return;
  }

  for (auto& parasite : *parasites_) {
    delete parasite;
  }
  parasites_->clear();
}

void SingleHostClonalParasitePopulations::add(ClonalParasitePopulation* blood_parasite) {
  blood_parasite->set_parasite_population(this);

  parasites_->push_back(blood_parasite);
  blood_parasite->set_index(parasites_->size() - 1);
  assert(parasites_->at(blood_parasite->index()) == blood_parasite);
}

void SingleHostClonalParasitePopulations::remove(ClonalParasitePopulation* blood_parasite) {
  remove(blood_parasite->index());
}

void SingleHostClonalParasitePopulations::remove(const int& index) {
  ClonalParasitePopulation* bp = parasites_->at(index);
  //    std::cout << parasites_.size() << std::endl;
  if (bp->index() != index) {
    LOG(FATAL) << "Incorrect index when remove parasite from SingleHostClonalParasitePopulations: " << bp->index()
               << "-" << index << "-" << parasites_->at(index)->index() << std::endl;
    assert(bp->index() == index);
  }

  parasites_->back()->set_index(index);
  parasites_->at(index) = parasites_->back();
  parasites_->pop_back();
  bp->set_index(-1);

  bp->set_parasite_population(nullptr);

  delete bp;
  bp = nullptr;
}

// return normalized parasite density and log 10 total infectious density
void SingleHostClonalParasitePopulations::get_parasites_profiles(std::vector<double>& relative_parasite_density,
                                                                 double& log10_total_relative_density) const {
  auto i = 0;

  while ((i < parasites_->size())
         && (NumberHelpers::is_equal(parasites_->at(i)->get_log10_relative_density(),
                                     ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY))) {
    relative_parasite_density[i] = 0.0;
    i++;
  }

  if (i == parasites_->size()) {
    log10_total_relative_density = ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
    return;
  }

  log10_total_relative_density = parasites_->at(i)->get_log10_relative_density();
  relative_parasite_density[i] = (log10_total_relative_density);

  for (auto j = i + 1; j < parasites_->size(); j++) {
    const auto log10_relative_density = parasites_->at(j)->get_log10_relative_density();

    if (NumberHelpers::is_enot_qual(log10_relative_density, ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY)) {
      relative_parasite_density[j] = (log10_relative_density);
      log10_total_relative_density += log10(pow(10, log10_relative_density - log10_total_relative_density) + 1);
    } else {
      relative_parasite_density[j] = (0.0);
    }
  }

  // normalize
  for (auto j = 0; j < parasites_->size(); j++) {
    if (NumberHelpers::is_enot_qual(relative_parasite_density[j], 0.0)) {
      relative_parasite_density[j] = pow(10, relative_parasite_density[j] - log10_total_relative_density);
    }
  }
}

double SingleHostClonalParasitePopulations::get_log10_total_relative_density() {
  auto i = 0;

  while ((i < parasites_->size())
         && (NumberHelpers::is_equal(parasites_->at(i)->get_log10_relative_density(),
                                     ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY))) {
    i++;
  }

  if (i == parasites_->size()) {
    return ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
  }

  auto log10_total_relative_density = parasites_->at(i)->get_log10_relative_density();

  for (int j = i + 1; j < parasites_->size(); j++) {
    const auto log10_relative_density = parasites_->at(j)->get_log10_relative_density();

    if (NumberHelpers::is_enot_qual(log10_relative_density, ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY)) {
      log10_total_relative_density += log10(pow(10, log10_relative_density - log10_total_relative_density) + 1);
    }
  }

  return log10_total_relative_density;
}

int SingleHostClonalParasitePopulations::latest_update_time() const {
  return person_->latest_update_time();
}

int SingleHostClonalParasitePopulations::size() {
  return parasites_->size();
}

bool SingleHostClonalParasitePopulations::contain(ClonalParasitePopulation* blood_parasite) {
  for (auto& parasite : *parasites_) {
    if (blood_parasite == parasite) {
      return true;
    }
  }

  return false;
}

void SingleHostClonalParasitePopulations::change_all_parasite_update_function(ParasiteDensityUpdateFunction* from,
                                                                              ParasiteDensityUpdateFunction* to) const {
  for (auto* parasite : *parasites_) {
    if (parasite->update_function() == from) {
      parasite->set_update_function(to);
    }
  }
}

void SingleHostClonalParasitePopulations::update() const {
  for (auto* bp : *parasites_) {
    bp->update();
  }
  //    std::vector<BloodParasite*>(*parasites_).swap(*parasites_);
}

void SingleHostClonalParasitePopulations::clear_cured_parasites() {
  //    std::vector<int> cured_parasites_index;
  for (int i = parasites_->size() - 1; i >= 0; i--) {
    if (parasites_->at(i)->last_update_log10_parasite_density()
        <= Model::CONFIG->parasite_density_level().log_parasite_density_cured + 0.00001) {
      remove(i);
    }
  }
}

void SingleHostClonalParasitePopulations::update_by_drugs(DrugsInBlood* drugs_in_blood) const {
  for (auto& blood_parasite : *parasites_) {
    auto* new_genotype = blood_parasite->genotype();

    double percent_parasite_remove = 0;
    for (auto& [drug_id, drug] : *drugs_in_blood->drugs()) {
      const auto p = Model::RANDOM->random_flat(0.0, 1.0);

      if (p < drug->get_mutation_probability()) {
        // select all locus
        // TODO: rework here
        // remember to use mask to turn on and off mutation location
        // for a specific time
        Genotype* candidate_genotype =
            new_genotype->perform_mutation_by_drug(Model::CONFIG, Model::RANDOM, drug->drug_type());

        if (candidate_genotype->get_EC50_power_n(drug->drug_type())
            > new_genotype->get_EC50_power_n(drug->drug_type())) {
          // higher EC50^n means lower efficacy then allow mutation occur
          new_genotype = candidate_genotype;
        }
        if (new_genotype != blood_parasite->genotype()) {
          // mutation occurs
          Model::DATA_COLLECTOR->record_1_mutation(person_->location(), blood_parasite->genotype(), new_genotype);
//          LOG(TRACE) << Model::SCHEDULER->current_time() << "\t" << blood_parasite->genotype()->genotype_id << "\t"
//                     << new_genotype->genotype_id << "\t"
//                     << blood_parasite->genotype()->get_EC50_power_n(drug->drug_type()) << "\t"
//                     << new_genotype->get_EC50_power_n(drug->drug_type());
          blood_parasite->set_genotype(new_genotype);
        }
      }

      const auto p_temp = drug->get_parasite_killing_rate(blood_parasite->genotype()->genotype_id);
      percent_parasite_remove = percent_parasite_remove + p_temp - percent_parasite_remove * p_temp;
    }
    if (percent_parasite_remove > 0) {
      blood_parasite->perform_drug_action(percent_parasite_remove);
    }
  }
}

bool SingleHostClonalParasitePopulations::has_detectable_parasite() const {
  for (auto& parasite : *parasites_) {
    if (parasite->last_update_log10_parasite_density()
        >= Model::CONFIG->parasite_density_level().log_parasite_density_detectable_pfpr) {
      return true;
    }
  }
  return false;
}

bool SingleHostClonalParasitePopulations::is_gametocytaemic() const {
  for (auto& parasite : *parasites_) {
    if (parasite->gametocyte_level() > 0) {
      return true;
    }
  }
  return false;
}
