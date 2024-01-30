/* 
 * File:   DrugsInBlood.cpp
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
  
  person_->person_index_gpu->h_person_update_info()[person_->index()].drug_in_blood_size = drugs_->size();
  if(person_->person_index_gpu->h_person_update_info()[person_->index()].drug_in_blood_size > 0){
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_starting_value[typeID] = drugs_->at(typeID)->starting_value();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_dosing_days[typeID] = drugs_->at(typeID)->dosing_days();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_last_update_value[typeID] = drugs_->at(typeID)->last_update_value();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_last_update_time[typeID] = drugs_->at(typeID)->last_update_time();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_start_time[typeID] = drugs_->at(typeID)->start_time();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_end_time[typeID] = drugs_->at(typeID)->end_time();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_half_life[typeID] = drugs_->at(typeID)->drug_type()->drug_half_life();
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_rand_uniform_1[typeID] = 0.1;//Model::RANDOM->random_uniform_double(-0.2, 0.2);
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_rand_uniform_2[typeID] = 0.1;//Model::RANDOM->random_uniform_double(0, 0.1);
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_in_blood_type_id[person_->person_index_gpu->h_person_update_info()[person_->index()].drug_in_blood_size-1] = typeID;
    person_->person_index_gpu->h_person_update_info()[person_->index()].drug_epsilon = std::numeric_limits<double>::epsilon();
  }
  if(person_->index() >= 1040 && person_->index() <= 1045)
    printf("GPU::DrugsInBlood::add_drug() person_->index() = %d %d %d %d %d %f %f %f\n",
           person_->index(),drugs_->size(),typeID,
           person_->person_index_gpu->h_person_update_info()[person_->index()].drug_start_time[typeID],
           person_->person_index_gpu->h_person_update_info()[person_->index()].drug_last_update_time[typeID],
           person_->person_index_gpu->h_person_update_info()[person_->index()].drug_starting_value[typeID],
           person_->person_index_gpu->h_person_update_info()[person_->index()].drug_last_update_value[typeID],
           person_->person_index_gpu->h_person_update_info()[person_->index()].drug_half_life[typeID]);
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
    if(person_->index() >= 1040 && person_->index() <= 1045)
      printf("GPU::DrugsInBlood::clear_cut_off_drugs_by_event() before person_->index() = %d %d\n",
             person_->index(),drugs_->size());
    for (auto pos = drugs_->begin(); pos!=drugs_->end();) {
      //if (pos->second->lastUpdateValue <= 0.1) {
      //Cut off at 10%
      if (pos->second->last_update_value() <= DRUG_CUT_OFF_VALUE) {
        //if drug is astermisinin then deActive Gametocyte

        //                person->cancelClearDrugEvent(pos->first, eventID);
        person_->person_index_gpu->h_person_update_info()[person_->index()].drug_in_blood_type_id[pos->second->drug_type()->id()] = -1;
        delete pos->second;
        drugs_->erase(pos++);
      } else {
        ++pos;
      }
    }
    if(person_->index() >= 1040 && person_->index() <= 1045)
      printf("GPU::DrugsInBlood::clear_cut_off_drugs_by_event() after person_->index() = %d %d\n",
             person_->index(),drugs_->size());
  }
}
