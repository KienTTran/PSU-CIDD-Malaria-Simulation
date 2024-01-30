/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TherapyBuilder.cpp
 * Author: Merlin
 * 
 * Created on August 23, 2017, 2:43 PM
 */

#include "TherapyBuilder.cuh"
#include "Therapy.cuh"
#include "SCTherapy.cuh"
#include "MACTherapy.cuh"

GPU::TherapyBuilder::TherapyBuilder() = default;

GPU::TherapyBuilder::~TherapyBuilder() = default;

GPU::Therapy* GPU::TherapyBuilder::build(const YAML::Node &ns, const int &t_id) {
  GPU::Therapy* t = nullptr;
  if (ns["drug_id"]) {
    t = new GPU::SCTherapy();
    t->set_id(t_id);
    for (auto i = 0; i < ns["drug_id"].size(); i++) {
      const auto drug_id = ns["drug_id"][i].as<int>();
      //        std::cout << therapy_id << "-" << drug_id << std::endl;
      dynamic_cast<GPU::SCTherapy*>(t)->add_drug(drug_id);
    }
    for (auto i = 0; i < ns["dosing_days"].size(); i++) {
      const auto dosing_days = ns["dosing_days"][i].as<int>();
      dynamic_cast<GPU::SCTherapy*>(t)->dosing_day.push_back(dosing_days);
    }
  } else {
    if (ns["therapy_ids"]) {
      t = new GPU::MACTherapy();
      t->set_id(t_id);

      for (int i = 0; i < ns["therapy_ids"].size(); i++) {
        auto therapy_id = ns["therapy_ids"][i].as<int>();
        //        std::cout << therapy_id << "-" << drug_id << std::endl;
        dynamic_cast<GPU::MACTherapy*>(t)->add_therapy_id(therapy_id);
      }
      for (int i = 0; i < ns["regimen"].size(); i++) {
        auto starting_day = ns["regimen"][i].as<int>();
        //        std::cout << therapy_id << "-" << drug_id << std::endl;
        dynamic_cast<GPU::MACTherapy*>(t)->add_schedule(starting_day);
      }
    }
  }

  if (ns["name"]){
    t->set_name(ns["name"].as<std::string>());
  }

  return t;
}
