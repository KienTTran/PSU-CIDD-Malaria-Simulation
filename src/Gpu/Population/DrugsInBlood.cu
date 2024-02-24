/* 
 * File:   DrugsInBlood.cu
 * Author: Merlin
 * 
 * Created on July 31, 2013, 1:47 PM
 */

#include "DrugsInBlood.cuh"
#include "Gpu/Therapies/Drug.cuh"
#include "Gpu/Events/Event.cuh"
#include "Gpu/Therapies/DrugType.cuh"
#include "Person.cuh"
#include "Helpers/ObjectHelpers.h"
#include "Core/TypeDef.h"
#include "Core/Random.h"
#include "Model.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"

#ifndef DRUG_CUT_OFF_VALUE
#define DRUG_CUT_OFF_VALUE 0.1
#endif

GPU::DrugsInBlood::DrugsInBlood(GPU::Person *person) : person_(person), drugs_(nullptr) {}

void GPU::DrugsInBlood::init() {
  drugs_ = new GPUDrugPtrMap();
}

GPU::DrugsInBlood::~DrugsInBlood() {
  if (drugs_!=nullptr) {
    clear();
    ObjectHelpers::delete_pointer<GPUDrugPtrMap>(drugs_);
  }
}

GPU::Drug *GPU::DrugsInBlood::add_drug(GPU::Drug *drug) {
  int typeID = drug->drug_type()->id();
  if (!is_drug_in_blood(typeID)) {
    drug->set_person_drugs(this);
    drugs_->insert(std::pair<int, GPU::Drug *>(typeID, drug));
  } else {
    //already have it
    drugs_->at(typeID)->set_starting_value(drug->starting_value());
    drugs_->at(typeID)->set_dosing_days(drug->dosing_days());
    drugs_->at(typeID)->set_last_update_value(drug->last_update_value());
    drugs_->at(typeID)->set_last_update_time(drug->last_update_time());
    drugs_->at(typeID)->set_start_time(drug->start_time());
    drugs_->at(typeID)->set_end_time(drug->end_time());
    delete drug;
  }

  auto &person_update_info = person_->person_index_gpu->h_person_update_info()[person_->index()];
  person_update_info.add_drug_in_blood(drugs_,typeID);

  LOG_IF(person_->index() >= 1000 && person_->index() <= 1085,INFO)
    << fmt::format("GPU::DrugsInBlood::add_drug() person_->index() = {} {} {} {} {} {} {} {}",
           person_->index(),drugs_->size(),typeID,
           person_update_info.drug_start_time[typeID],
           person_update_info.drug_last_update_time[typeID],
           person_update_info.drug_starting_value[typeID],
           person_update_info.drug_last_update_value[typeID],
           person_update_info.drug_half_life[typeID]);

  return drugs_->at(typeID);
}

bool GPU::DrugsInBlood::is_drug_in_blood(GPU::DrugType *drug_type) const {
  return is_drug_in_blood(drug_type->id());
}

bool GPU::DrugsInBlood::is_drug_in_blood(const int drugTypeID) const {
  return drugs_->find(drugTypeID)!=drugs_->end();
}

void GPU::DrugsInBlood::remove_drug(GPU::Drug *drug) const {
  remove_drug(drug->drug_type()->id());
}

void GPU::DrugsInBlood::remove_drug(const int &drug_type_id) const {
  auto it = drugs_->find(drug_type_id);

  if (it==drugs_->end()) {
    return;
  }

  delete it->second;
  drugs_->erase(it);
}

GPU::Drug *GPU::DrugsInBlood::get_drug(const int &type_id) const {
  if (!is_drug_in_blood(type_id))
    return nullptr;

  return drugs_->at(type_id);
}

std::size_t GPU::DrugsInBlood::size() const {
  return drugs_->size();
}

void GPU::DrugsInBlood::clear() const {
  if (drugs_->empty()) return;

  for (auto &drug : *drugs_) {
    delete drug.second;
  }
  drugs_->clear();
}

void GPU::DrugsInBlood::update() const {
  for (auto &drug : *drugs_) {
    drug.second->update();
  }
}

void GPU::DrugsInBlood::clear_cut_off_drugs_by_event(GPU::Event *event) const {
  if (!drugs_->empty()) {
    LOG_IF(person_->index() >= 1000 && person_->index() <= 1085,INFO)
      << fmt::format("GPU::DrugsInBlood::clear_cut_off_drugs_by_event() before person_->index() = {} {}",
             person_->index(),drugs_->size());
    for (auto pos = drugs_->begin(); pos!=drugs_->end();) {
      //if (pos->second->lastUpdateValue <= 0.1) {
      //Cut off at 10%
      if (pos->second->last_update_value() <= DRUG_CUT_OFF_VALUE) {
        //if drug is astermisinin then deActive Gametocyte

        //                person->cancelClearDrugEvent(pos->first, eventID);
        delete pos->second;
        drugs_->erase(pos++);
      } else {
        ++pos;
      }
    }
    LOG_IF(person_->index() >= 1000 && person_->index() <= 1085,INFO)
      << fmt::format("GPU::DrugsInBlood::clear_cut_off_drugs_by_event() after person_->index() = {} {}",
             person_->index(),drugs_->size());
  }
}
