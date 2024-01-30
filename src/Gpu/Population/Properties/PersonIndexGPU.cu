//
// Created by kient on 12/9/2023.
//

#include "Model.h"
#include "Core/Config/Config.h"
#include "PersonIndexGPU.cuh"
#include "Gpu/Population/SingleHostClonalParasitePopulations.cuh"
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

//void GPU::PersonIndexGPU::add(GPU::Person *p,const int &location, const GPU::Person::HostStates &host_state, const int &age_class) {
////    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
////        << fmt::format("PersonIndexGPU add before loc {}->{} hs {}->{} ac {}->{} size {} index {}",
////                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
////                       h_persons().size(),
////                       p->GPU::PersonIndexGPUHandler::index());
//    /*
//     * This is re-add person so person attributes are set
//     * */
//    h_persons_.push_back(p);
//    p->set_index(h_persons_.size() - 1);
//    p->GPU::PersonIndexGPUHandler::set_index(h_persons().size() - 1);
//    p->set_id(h_persons_.size());
//    h_person_ids_.push_back(p->id());
////    printf("GPU::PersonIndexGPU::add2 loc %d hs %d ac %d person_index %d %d index2 %d id %d id2 %d\n",
////           location,host_state,age_class,
////           p->index(),
////           p->GPU::PersonIndexGPUHandler::index(),
////           p->person_update_info.person_index,
////           p->id(),
////           p->person_update_info.person_id);
//    h_person_host_states_.push_back(static_cast<int>(host_state));
//    h_person_ages_.push_back(p->age());
//    h_person_age_classes_.push_back(age_class);
//    h_person_locations_.push_back(location);
//    h_person_residence_locations_.push_back(p->residence_location());
//    if(Model::CONFIG->render_config().display_gui) {
//        if (p->residence_location() != location) {
//            p->generate_render_entity(location, true);
//        } else {
//            p->generate_render_entity(location);
//        }
//        h_person_models_[p->GPU::PersonIndexGPUHandler::index()] = p->model();
//        h_person_colors_[p->GPU::PersonIndexGPUHandler::index()] = p->color();
//    }
//    h_person_update_info_.push_back(p->person_update_info());
////    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
////        << fmt::format("PersonIndexGPU add after loc {}->{} hs {}->{} ac {}->{} size {} index {}",
////                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
////                       h_persons().size(),
////                       p->GPU::PersonIndexGPUHandler::index());
//}

GPU::PersonUpdateInfo GPU::PersonIndexGPU::init_person_update_info(GPU::Person *p) {
    GPU::PersonUpdateInfo person_update_info;
    person_update_info.person_index = p->GPU::PersonIndexGPUHandler::index();
    person_update_info.person_id = p->id();
    person_update_info.person_id = -1;
    person_update_info.person_index = -1;
    person_update_info.person_latest_update_time = -1;
    person_update_info.person_latest_immune_value = -1.0;
    /* for parasite update */
    person_update_info.parasite_size = -1;
    person_update_info.parasite_id = -1;
    person_update_info.parasite_current_index = -1;
    for(int i = 0; i < MAX_PARASITE_PER_PERSON; i++){
        person_update_info.parasite_update_function_type[i] = 0;
        person_update_info.parasite_last_update_log10_parasite_density[i] = GPU::ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;
        person_update_info.parasite_genotype_fitness_multiple_infection[i] = 1.0;
        person_update_info.parasite_gametocyte_level[i] = 0.0;
    }
    /* for drug update */
    for(int i = 0; i < MAX_DRUG_PER_ACT; i++){
        person_update_info.drug_in_blood_type_id[i] = -1;
    }
    return person_update_info;
}

void GPU::PersonIndexGPU::add(GPU::Person *p) {
    h_persons_.push_back(p);
    p->GPU::PersonIndexGPUHandler::set_index(h_persons_.size() - 1);
    p->set_index(p->GPU::PersonIndexGPUHandler::index());
    h_person_ids_.push_back(p->id());
    h_person_host_states_.push_back(static_cast<int>(p->host_state()));
    h_person_ages_.push_back(p->age());
    h_person_age_classes_.push_back(p->age_class());
    h_person_locations_.push_back(p->location());
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

    h_person_ids_[p->GPU::PersonIndexGPUHandler::index()] = h_person_ids_.back();
    h_person_ids_.pop_back();

    h_person_host_states_[p->GPU::PersonIndexGPUHandler::index()] = h_person_host_states_.back();
    h_person_host_states_.pop_back();

    h_person_ages_[p->GPU::PersonIndexGPUHandler::index()] = h_person_ages_.back();
    h_person_ages_.pop_back();

    h_person_age_classes_[p->GPU::PersonIndexGPUHandler::index()] = h_person_age_classes_.back();
    h_person_age_classes_.pop_back();

    h_person_locations_[p->GPU::PersonIndexGPUHandler::index()] = h_person_locations_.back();
    h_person_locations_.pop_back();

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
    switch (property) {
        case Person::LOCATION: {
//            printf("GPU::PersonIndexGPU::notify_change LOCATION\n");
//            change_property(p, *(int *) newValue, p->host_state(), p->age_class());
            p->location_ = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->location_ = *(int *) newValue;
            h_person_locations_[p->GPU::PersonIndexGPUHandler::index()] = *(int *) newValue;
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
//            change_property(p, p->location(), *(GPU::Person::HostStates *) newValue, p->age_class());
            p->host_state_ = *(GPU::Person::HostStates *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->host_state_ = *(GPU::Person::HostStates *) newValue;
            h_person_host_states_[p->GPU::PersonIndexGPUHandler::index()] = *(GPU::Person::HostStates *) newValue;
//            printf("GPU::PersonIndexGPU::notify_change HOST_STATE AFTER %d -> %d\n",p->host_state(),*(GPU::Person::HostStates *) newValue);
            break;
        }
        case Person::AGE: {
//            printf("GPU::PersonIndexGPU::notify_change AGE\n");
//            change_property(p, p->location(), p->host_state(), *(int *) newValue);
            p->age_ = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->age_ = *(int *) newValue;
            h_person_ages_[p->GPU::PersonIndexGPUHandler::index()] = *(int *) newValue;
            break;
        }
        case Person::AGE_CLASS: {
//            printf("GPU::PersonIndexGPU::notify_change AGE_CLASS\n");
//            change_property(p, p->location(), p->host_state(), *(int *) newValue);
            p->age_class_ = *(int *) newValue;
            h_persons_[p->GPU::PersonIndexGPUHandler::index()]->age_class_ = *(int *) newValue;
            h_person_age_classes_[p->GPU::PersonIndexGPUHandler::index()] = *(int *) newValue;
            break;
        }
        default:break;
    }
}

void GPU::PersonIndexGPU::update() {
//    printf("GPU::PersonIndexGPU::update() is called\n");
    h_persons_.shrink_to_fit();
    h_person_ids_.shrink_to_fit();
    h_person_host_states_.shrink_to_fit();
    h_person_ages_.shrink_to_fit();
    h_person_age_classes_.shrink_to_fit();
    h_person_locations_.shrink_to_fit();
    h_person_residence_locations_.shrink_to_fit();
    if(Model::CONFIG->render_config().display_gui) {
        h_person_models_.shrink_to_fit();
        h_person_colors_.shrink_to_fit();
    }
    h_person_update_info_.shrink_to_fit();
}
