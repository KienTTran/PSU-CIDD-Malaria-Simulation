/*
 * File:   ParasitePopulation.cpp
 * Author: Merlin
 *
 * Created on July 11, 2013, 1:53 PM
 */

#include "SingleHostClonalParasitePopulations.cuh"

#include <cmath>

#include "ClonalParasitePopulation.cuh"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Helpers/NumberHelpers.h"
#include "Helpers/ObjectHelpers.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Model.h"
#include "Therapies/Drug.h"
#include "Gpu/Parasites/Genotype.cuh"
#include "DrugsInBlood.cuh"
#include "Person.cuh"

__device__ __host__ GPU::SingleHostClonalParasitePopulations::SingleHostClonalParasitePopulations(GPU::Person* person)
    : person_(person),
      parasites_(nullptr),
      log10_total_infectious_denstiy(-1000.0) {}

void GPU::SingleHostClonalParasitePopulations::init() {
  parasites_ = new std::vector<GPU::ClonalParasitePopulation*>();
}

GPU::SingleHostClonalParasitePopulations::~SingleHostClonalParasitePopulations() {
  if (parasites_ != nullptr) {
    clear();
    ObjectHelpers::delete_pointer<std::vector<GPU::ClonalParasitePopulation*>>(parasites_);
  }

  person_ = nullptr;
}

void GPU::SingleHostClonalParasitePopulations::clear() {
  if (parasites_->empty()) {
    return;
  }

  for (auto& parasite : *parasites_) {
    delete parasite;
  }
  parasites_->clear();
}

void GPU::SingleHostClonalParasitePopulations::add(GPU::ClonalParasitePopulation* blood_parasite) {
  blood_parasite->set_parasite_population(this);
  parasites_->push_back(blood_parasite);
  blood_parasite->set_index(parasites_->size() - 1);
  assert(parasites_->at(blood_parasite->index()) == blood_parasite);
}

void GPU::SingleHostClonalParasitePopulations::remove(GPU::ClonalParasitePopulation* blood_parasite) {
  remove(blood_parasite->index());
}

void GPU::SingleHostClonalParasitePopulations::remove(const int& index) {
    GPU::ClonalParasitePopulation* bp = parasites_->at(index);
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

int GPU::SingleHostClonalParasitePopulations::latest_update_time() const {
  return person_->latest_update_time();
}

int GPU::SingleHostClonalParasitePopulations::size() {
  return parasites_->size();
}

bool GPU::SingleHostClonalParasitePopulations::contain(GPU::ClonalParasitePopulation* blood_parasite) {
  for (auto& parasite : *parasites_) {
    if (blood_parasite == parasite) {
      return true;
    }
  }

  return false;
}

void GPU::SingleHostClonalParasitePopulations::change_all_parasite_update_function(GPU::ParasiteDensityUpdateFunction* from,
                                                                                   GPU::ParasiteDensityUpdateFunction* to) const {
  for (auto* parasite : *parasites_) {
    if (parasite->update_function() == from) {
      parasite->set_update_function(to);
      parasite->set_gpu_update_function(to);
    }
  }
}

void GPU::SingleHostClonalParasitePopulations::update() {
  for (auto* bp : *parasites_) {
    bp->update();
  }
}

void GPU::SingleHostClonalParasitePopulations::clear_cured_parasites() {
  log10_total_infectious_denstiy = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
  //    std::vector<int> cured_parasites_index;
  for (int i = parasites_->size() - 1; i >= 0; i--) {
    if (parasites_->at(i)->last_update_log10_parasite_density()
        <= Model::CONFIG->parasite_density_level().log_parasite_density_cured + 0.00001) {
      remove(i);
    } else {
      if (log10_total_infectious_denstiy == GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY) {
        log10_total_infectious_denstiy = (*parasites_)[i]->get_log10_infectious_density();
      } else {
        log10_total_infectious_denstiy +=
            log10(pow(10, (*parasites_)[i]->get_log10_infectious_density() - log10_total_infectious_denstiy) + 1);
      }
    }
  }
}

void GPU::SingleHostClonalParasitePopulations::update_by_drugs(GPU::DrugsInBlood* drugs_in_blood) const {
  for (auto& blood_parasite : *parasites_) {
    auto* new_genotype = blood_parasite->genotype();

    double percent_parasite_remove = 0;
    for (auto& [drug_id, drug] : *drugs_in_blood->drugs()) {
      // select all locus
      // remember to use mask to turn on and off mutation location
      // for a specific time
      GPU::Genotype* candidate_genotype = new_genotype->perform_mutation_by_drug(Model::CONFIG, Model::RANDOM, drug->drug_type(),
                                                                            Model::CONFIG->mutation_probability_by_locus());

      if (candidate_genotype->get_EC50_power_n(drug->drug_type())
          > new_genotype->get_EC50_power_n(drug->drug_type())) {
        // higher EC50^n means lower efficacy then allow mutation occur
        new_genotype = candidate_genotype;
      }
      if (new_genotype != blood_parasite->genotype()) {
        // mutation occurs
        Model::GPU_DATA_COLLECTOR->record_1_mutation(person_->location(), blood_parasite->genotype(), new_genotype);
        //          LOG(TRACE) << Model::GPU_SCHEDULER->current_time() << "\t" << blood_parasite->genotype()->genotype_id
        //          << "\t"
        //                     << new_genotype->genotype_id << "\t"
        //                     << blood_parasite->genotype()->get_EC50_power_n(drug->drug_type()) << "\t"
        //                     << new_genotype->get_EC50_power_n(drug->drug_type());
        blood_parasite->set_genotype(new_genotype);
      }

      const auto p_temp = drug->get_parasite_killing_rate(blood_parasite->genotype()->genotype_id);
      percent_parasite_remove = percent_parasite_remove + p_temp - percent_parasite_remove * p_temp;
    }
    if (percent_parasite_remove > 0) {
      blood_parasite->perform_drug_action(percent_parasite_remove);
    }
  }
}

bool GPU::SingleHostClonalParasitePopulations::has_detectable_parasite() const {
  for (auto& parasite : *parasites_) {
    if (parasite->last_update_log10_parasite_density()
        >= Model::CONFIG->parasite_density_level().log_parasite_density_detectable_pfpr) {
      return true;
    }
  }
  return false;
}

bool GPU::SingleHostClonalParasitePopulations::is_gametocytaemic() const {
  for (auto& parasite : *parasites_) {
    if (parasite->gametocyte_level() > 0) {
      return true;
    }
  }
  return false;
}
