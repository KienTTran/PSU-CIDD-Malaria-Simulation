/* 
 * File:   DrugDatabase.cpp
 * Author: nguyentran
 * 
 * Created on June 3, 2013, 11:05 AM
 */


#include "DrugDatabase.cuh"

GPU::DrugDatabase::DrugDatabase() = default;

GPU::DrugDatabase::~DrugDatabase() {

  for (auto &i : *this) {
    delete i.second;
  }
  this->clear();
}

void GPU::DrugDatabase::add(GPU::DrugType *dt) {
  (*this)[dt->id()] = dt;
}

