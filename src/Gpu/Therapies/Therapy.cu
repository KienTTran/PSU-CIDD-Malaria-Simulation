/* 
 * File:   Therapy.cpp
 * Author: nguyentran
 * 
 * Created on June 3, 2013, 7:50 PM
 */

#include "Therapy.cuh"

GPU::Therapy::Therapy() : id_{-1}, testing_day_{-1}, drug_ids{}, name_{""} {
}

GPU::Therapy::~Therapy() = default;

void GPU::Therapy::add_drug(int drug_id) {
  drug_ids.push_back(drug_id);
}

//int GPU::Therapy::get_therapy_duration(int dosing_day) {
//    int result = 0;
//
//    for (int i = 0; i < drug_ids_.size(); i++) {
//        DrugType* dt = Model::CONFIG->gpu_drug_db()->get(drug_ids_[i]);
//        if (!dt->is_artemisinin()) {
//            result = std::max<int>(dt->get_total_duration_of_drug_activity(dosing_day), result);
//        }
//    }
//    return result;
//}