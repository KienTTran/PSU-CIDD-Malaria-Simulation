//
// Created by kient on 12/9/2023.
//

#include "Model.h"
#include "Core/Config/Config.h"
#include "PersonIndexGPU.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
#include "Gpu/Population/ImmuneSystem.cuh"
#include "Gpu/Population/ImmuneComponent.cuh"
#include "Helpers/UniqueId.hxx"

GPU::PersonIndexGPU::PersonIndexGPU() {
    //Allocate memory to render from beginning to end of simulation because OGL need to be pre-allocated in order to avoid init data everytime
    if(Model::CONFIG->render_config().display_gui){
        h_person_models_ = TVector<glm::mat4>(Model::CONFIG->n_people_init()*Model::CONFIG->gpu_config().pre_allocated_mem_ratio);
        h_person_colors_ = TVector<glm::vec4>(Model::CONFIG->n_people_init()*Model::CONFIG->gpu_config().pre_allocated_mem_ratio);
    }
}

GPU::PersonIndexGPU::~PersonIndexGPU() {
}

GPU::PersonUpdateInfo GPU::PersonIndexGPU::init_person_update_info(GPU::Person *p) {
    GPU::PersonUpdateInfo person_update_info;
    person_update_info.person_index = p->index();
    person_update_info.person_id = p->id();
    person_update_info.person_age = p->age();
    person_update_info.person_age_class = p->age_class();
    person_update_info.person_host_state = static_cast<int>(p->host_state());
    person_update_info.person_location = p->location();
    person_update_info.person_residence_location = p->residence_location();
    person_update_info.person_birthday = p->birthday();
    person_update_info.person_moving_level = p->moving_level();
    person_update_info.person_innate_relative_biting_rate = p->innate_relative_biting_rate;
    person_update_info.person_current_relative_biting_rate = p->current_relative_biting_rate;
    person_update_info.person_using_age_dependent_biting_level = Model::CONFIG->using_age_dependent_bitting_level();
    person_update_info.person_latest_update_time = p->latest_update_time();
    person_update_info.person_has_drug_in_blood = false;
    /* for parasite update */
    person_update_info.limit_epsilon = std::numeric_limits<double>::epsilon(); //For device use
    person_update_info.LOG_ZERO_PARASITE_DENSITY = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
    person_update_info.parasites_size = -1;
    person_update_info.parasites_current_index = -1;
    for(int i = 0; i < MAX_PARASITE_PER_PERSON; i++){
        person_update_info.parasite_id[i] = -1;
        person_update_info.parasite_update_function_type[i] = 0;
        person_update_info.parasite_last_update_log10_parasite_density[i] = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
        person_update_info.parasite_genotype_fitness_multiple_infection[i] = 1.0;
        person_update_info.parasite_gametocyte_level[i] = 0.0;
        person_update_info.parasite_log10_infectious_density[i] = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
    }
    person_update_info.parasites_log10_total_infectious_density = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
    person_update_info.parasites_genotype_mutated_number = 0;
    /* for drug update */
    for(int i = 0; i < MAX_DRUG_PER_ACT; i++){
        person_update_info.drug_in_blood_type_id[i] = -1;
    }
    /* for immune system update*/
    person_update_info.immune_system_component_latest_value = p->immune_system()->get_lastest_immune_value();
    person_update_info.immune_system_component_type = p->immune_system()->immune_component()->type();
    person_update_info.immune_system_is_increased = p->immune_system()->increase();
    return person_update_info;
}

void GPU::PersonIndexGPU::update_person(){
    for(int i = 0; i < h_persons_.size(); i++){
        h_persons_[i]->set_host_state(static_cast<GPU::Person::HostStates>(h_person_update_info_[i].person_host_state));
        h_persons_[i]->immune_system()->set_increase(h_person_update_info_[i].immune_system_is_increased);
    }
}


void GPU::PersonIndexGPU::add(GPU::Person *p) {
    h_persons_.push_back(p);
    p->GPU::PersonIndexGPUHandler::set_index(h_persons_.size() - 1);
    p->set_index(p->GPU::PersonIndexGPUHandler::index());
    h_person_residence_locations_.push_back(p->residence_location());
    if(Model::CONFIG->render_config().display_gui) {
        p->generate_render_entity(p->location());
        h_person_models_[p->GPU::PersonIndexGPUHandler::index()] = p->model();
        h_person_colors_[p->GPU::PersonIndexGPUHandler::index()] = p->color();
    }
    h_person_update_info_.push_back(init_person_update_info(p));
}

void GPU::PersonIndexGPU::remove(GPU::Person *p) {

    remove_without_set_index(p);
    p->GPU::PersonIndexGPUHandler::set_index(-1);
    p->set_index(-1);
    //    delete p;
    //    p = nullptr;
}
/*
 *
 * Set the person to remove to last person
 * Remove the last person
 * Insert the property changed person to last index
 * */
void GPU::PersonIndexGPU::remove_without_set_index(GPU::Person *p) {
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU remove before loc {} hs {} ac {} size {} size2 {} index {} id {}",
//                             p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                             h_person_host_states_.size(),h_persons_.size(),
//                             p->GPU::PersonIndexGPUHandler::index(),p->id());

    h_persons_.back()->GPU::PersonIndexGPUHandler::set_index(p->GPU::PersonIndexGPUHandler::index());
    h_persons_.back()->set_index(p->GPU::PersonIndexGPUHandler::index());
    h_persons_.back()->set_id(p->id());

    h_person_residence_locations_[p->GPU::PersonIndexGPUHandler::index()] = h_person_residence_locations_.back();
    h_person_residence_locations_.pop_back();

    if(Model::CONFIG->render_config().display_gui){
        h_person_models_[p->GPU::PersonIndexGPUHandler::index()] = h_person_models_.back();
//    h_person_models_.pop_back(); //No remove because of OGL buffer need to be set to cover all populations for all time

        h_person_colors_[p->GPU::PersonIndexGPUHandler::index()] = h_person_colors_.back();
//    h_person_colors_.pop_back(); //No remove because of OGL buffer need to be set to cover all populations for all time
    }

    h_person_update_info_[p->GPU::PersonIndexGPUHandler::index()] = h_person_update_info_.back();
    h_person_update_info_.pop_back();

    //move the last element to current position and remove the last holder
    h_persons_[p->GPU::PersonIndexGPUHandler::index()] = h_persons_.back();
    h_persons_.pop_back();
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU remove after loc {} hs {} ac {} size {} size2 {} index {} id {}",
//                       p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                       h_person_host_states_.size(),h_persons_.size(),
//                       p->GPU::PersonIndexGPUHandler::index(),p->id());
}

std::size_t GPU::PersonIndexGPU::size() const {
    return h_persons_.size();
}

/*
 * Here we only update the person index,
 * no need to remove and add like N-dim index
 * */
void GPU::PersonIndexGPU::notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue,
                                  const void *newValue) {
    /*
     * Here using p->[property]_ to set value instead of using set_[property]_
     * otherwise it will be looped
     * */
    switch (property) {
        case Person::LOCATION: {
//            printf("GPU::PersonIndexGPU::notify_change LOCATION\n");
            p->location_ = *(int *) newValue;
            p->person_index_gpu->h_person_update_info_[p->index_].person_location = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->location_ = *(int *) newValue;
            if(Model::CONFIG->render_config().display_gui) {
                if (p->residence_location() != *(int *) newValue) {
                    p->generate_render_entity(*(int *) newValue, true);
                } else {
                    p->generate_render_entity(*(int *) newValue);
                }
                p->generate_render_entity(p->location());
                h_person_models_[p->GPU::PersonIndexGPUHandler::index()] = p->model();
                h_person_colors_[p->GPU::PersonIndexGPUHandler::index()] = p->color();
            }
            break;
        }
        case Person::HOST_STATE: {
//            printf("GPU::PersonIndexGPU::notify_change HOST_STATE BEFORE %d -> %d\n",p->host_state(),*(GPU::Person::HostStates *) newValue);
            p->host_state_ = *(GPU::Person::HostStates *) newValue;
            p->person_index_gpu->h_person_update_info_[p->index_].person_host_state = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->host_state_ = *(GPU::Person::HostStates *) newValue;
//            printf("GPU::PersonIndexGPU::notify_change HOST_STATE AFTER %d -> %d\n",p->host_state(),*(GPU::Person::HostStates *) newValue);
            break;
        }
        case Person::AGE: {
//            printf("GPU::PersonIndexGPU::notify_change AGE\n");
            p->age_ = *(int *) newValue;
            p->person_index_gpu->h_person_update_info_[p->index_].person_age = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->age_ = *(int *) newValue;
            break;
        }
        case Person::AGE_CLASS: {
//            printf("GPU::PersonIndexGPU::notify_change AGE_CLASS\n");
            p->age_class_ = *(int *) newValue;
            p->person_index_gpu->h_person_update_info_[p->index_].person_age_class = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->age_class_ = *(int *) newValue;
            break;
        }
        case Person::MOVING_LEVEL: {
//            printf("GPU::PersonIndexGPU::notify_change MOVING_LEVEL\n");
            p->moving_level_ = *(int *) newValue;
            p->person_index_gpu->h_person_update_info_[p->index_].person_moving_level = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->moving_level_ = *(int *) newValue;
            break;
        }
        default:break;
    }
}

void GPU::PersonIndexGPU::update() {
//    printf("GPU::PersonIndexGPU::update() is called\n");
    h_persons_.shrink_to_fit();
//    h_person_ids_.shrink_to_fit();
//    h_person_host_states_.shrink_to_fit();
//    h_person_ages_.shrink_to_fit();
//    h_person_age_classes_.shrink_to_fit();
//    h_person_locations_.shrink_to_fit();
    h_person_residence_locations_.shrink_to_fit();
    if(Model::CONFIG->render_config().display_gui) {
        h_person_models_.shrink_to_fit();
        h_person_colors_.shrink_to_fit();
    }
    h_person_update_info_.shrink_to_fit();
}
