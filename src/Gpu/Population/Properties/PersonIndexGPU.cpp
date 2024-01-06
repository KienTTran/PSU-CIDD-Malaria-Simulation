//
// Created by kient on 12/9/2023.
//

#include "Model.h"
#include "Core/Config/Config.h"
#include "PersonIndexGPU.h"

PersonIndexGPU::PersonIndexGPU() {
    //Allocate memory to render from beginning to end of simulation because OGL need to be pre-allocated in order to avoid init data everytime
    if(Model::CONFIG->render_config().display_gui){
        h_person_models_ = ThrustTVectorHost<glm::mat4>(Model::CONFIG->n_people_init()*Model::CONFIG->gpu_config().pre_allocated_mem_ratio);
        h_person_colors_ = ThrustTVectorHost<glm::vec4>(Model::CONFIG->n_people_init()*Model::CONFIG->gpu_config().pre_allocated_mem_ratio);
    }
}

PersonIndexGPU::~PersonIndexGPU() {
}

void PersonIndexGPU::add(Person *p,const int &location, const Person::HostStates &host_state, const int &age_class) {
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU add before loc {}->{} hs {}->{} ac {}->{} size {} index {}",
//                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                       h_persons().size(),
//                       p->PersonIndexGPUHandler::index());
    h_persons_.push_back(p);
    p->PersonIndexGPUHandler::set_index(h_persons_.size() - 1);
    h_person_ids_.push_back(p->id());
    h_person_host_states_.push_back(static_cast<int>(host_state));
    h_person_ages_.push_back(p->age());
    h_person_age_classes_.push_back(age_class);
    h_person_locations_.push_back(location);
    h_person_residence_locations_.push_back(p->residence_location());
    h_person_innate_relative_biting_rates_.push_back(p->innate_relative_biting_rate);
    h_person_current_relative_biting_rates_.push_back(p->current_relative_biting_rate);
    if(Model::CONFIG->render_config().display_gui) {
        if (p->residence_location() != location) {
            p->generate_render_entity(location, true);
        } else {
            p->generate_render_entity(location);
        }
        h_person_models_[p->PersonIndexGPUHandler::index()] = p->model();
        h_person_colors_[p->PersonIndexGPUHandler::index()] = p->color();
    }
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU add after loc {}->{} hs {}->{} ac {}->{} size {} index {}",
//                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                       h_persons().size(),
//                       p->PersonIndexGPUHandler::index());
}

void PersonIndexGPU::add(Person *p) {
    h_persons_.push_back(p);
    p->PersonIndexGPUHandler::set_index(h_persons_.size() - 1);
    h_person_ids_.push_back(p->id());
    h_person_host_states_.push_back(static_cast<int>(p->host_state()));
    h_person_ages_.push_back(p->age());
    h_person_age_classes_.push_back(p->age_class());
    h_person_locations_.push_back(p->location());
    h_person_residence_locations_.push_back(p->residence_location());
    h_person_innate_relative_biting_rates_.push_back(p->innate_relative_biting_rate);
    h_person_current_relative_biting_rates_.push_back(p->current_relative_biting_rate);
    if(Model::CONFIG->render_config().display_gui) {
        p->generate_render_entity(p->location());
        h_person_models_[p->PersonIndexGPUHandler::index()] = p->model();
        h_person_colors_[p->PersonIndexGPUHandler::index()] = p->color();
    }
}

void PersonIndexGPU::remove(Person *p) {

    remove_without_set_index(p);
    p->PersonIndexGPUHandler::set_index(-1);
    //    delete p;
    //    p = nullptr;
}

void PersonIndexGPU::remove_without_set_index(Person *p) {
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU remove before loc {} hs {} ac {} size {} size2 {} index {} id {}",
//                             p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                             h_person_host_states_.size(),h_persons_.size(),
//                             p->PersonIndexGPUHandler::index(),p->id());
    h_persons_.back()->PersonIndexGPUHandler::set_index(p->PersonIndexGPUHandler::index());
    h_person_ids_[p->PersonIndexGPUHandler::index()] = h_person_ids_[h_persons_.size() - 1];
    h_person_ids_.pop_back();

    h_person_host_states_[p->PersonIndexGPUHandler::index()] = h_person_host_states_[h_persons_.size() - 1];
    h_person_host_states_.pop_back();

    h_person_ages_[p->PersonIndexGPUHandler::index()] = h_person_ages_[h_persons_.size() - 1];
    h_person_ages_.pop_back();

    h_person_age_classes_[p->PersonIndexGPUHandler::index()] = h_person_age_classes_[h_persons_.size() - 1];
    h_person_age_classes_.pop_back();

    h_person_locations_[p->PersonIndexGPUHandler::index()] = h_person_locations_[h_persons_.size() - 1];
    h_person_locations_.pop_back();

    h_person_residence_locations_[p->PersonIndexGPUHandler::index()] = h_person_residence_locations_[h_persons_.size() - 1];
    h_person_residence_locations_.pop_back();

    h_person_innate_relative_biting_rates_[p->PersonIndexGPUHandler::index()] = h_person_innate_relative_biting_rates_[h_persons_.size() - 1];
    h_person_innate_relative_biting_rates_.pop_back();

    h_person_current_relative_biting_rates_[p->PersonIndexGPUHandler::index()] = h_person_current_relative_biting_rates_[h_persons_.size() - 1];
    h_person_current_relative_biting_rates_.pop_back();

    if(Model::CONFIG->render_config().display_gui){
        h_person_models_[p->PersonIndexGPUHandler::index()] = h_person_models_[h_persons_.size() - 1];
//    h_person_models_.pop_back(); //No remove because of OGL buffer need to be set to cover all populations for all time

        h_person_colors_[p->PersonIndexGPUHandler::index()] = h_person_colors_[h_persons_.size() - 1];
//    h_person_colors_.pop_back(); //No remove because of OGL buffer need to be set to cover all populations for all time
    }

    //move the last element to current position and remove the last holder
    h_persons_[p->PersonIndexGPUHandler::index()] = h_persons_.back();
    h_persons_.pop_back();
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU remove after loc {} hs {} ac {} size {} size2 {} index {} id {}",
//                       p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                       h_person_host_states_.size(),h_persons_.size(),
//                       p->PersonIndexGPUHandler::index(),p->id());
}

std::size_t PersonIndexGPU::size() const {
    return h_persons_.size();
}

void PersonIndexGPU::notify_change(Person *p, const Person::Property &property, const void *oldValue,
                                  const void *newValue) {
    switch (property) {
        case Person::LOCATION: {
            change_property(p, *(int *) newValue, p->host_state(), p->age_class());
        }
            break;
        case Person::HOST_STATE: {
            change_property(p, p->location(), *(Person::HostStates *) newValue, p->age_class());
        }
            break;
        case Person::AGE_CLASS: {
            change_property(p, p->location(), p->host_state(), *(int *) newValue);
        }
            break;
        default:break;
    }
}

void PersonIndexGPU::change_property(Person *p, const int &location,const Person::HostStates &host_state, const int &age_class) {
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << "PersonIndexGPU Changing " << p->location() << "-" << static_cast<int>(p->host_state()) << "-" << p->age_class() << " to "
//        << location << "-" << static_cast<int>(host_state) << "-" << age_class;
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU change_property before remove-add loc {}->{} hs {}->{} ac {}->{} size {} index {}",
//                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                       h_persons().size(),
//                       p->PersonIndexGPUHandler::index());
    remove_without_set_index(p); //to save 1 set and improve performance since the index of p will changed when add
    //add to new position
    add(p,location,host_state,age_class);
//    LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//        << fmt::format("PersonIndexGPU change_property after remove-add loc {}->{} hs {}->{} ac {}->{} size {} index {}",
//                       p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                       h_persons().size(),
//                       p->PersonIndexGPUHandler::index());
}

void PersonIndexGPU::update() {
    printf("PersonIndexGPU::update() is called\n");
    h_persons_.shrink_to_fit();
    h_person_ids_.shrink_to_fit();
    h_person_host_states_.shrink_to_fit();
    h_person_ages_.shrink_to_fit();
    h_person_age_classes_.shrink_to_fit();
    h_person_locations_.shrink_to_fit();
    h_person_residence_locations_.shrink_to_fit();
    h_person_innate_relative_biting_rates_.shrink_to_fit();
    h_person_current_relative_biting_rates_.shrink_to_fit();
    if(Model::CONFIG->render_config().display_gui) {
        h_person_models_.shrink_to_fit();
        h_person_colors_.shrink_to_fit();
    }
}
