//
// Created by kient on 12/9/2023.
//

#ifndef MASS_PERSONINDEXGPU_H
#define MASS_PERSONINDEXGPU_H

#include "Population/Properties/PersonIndex.h"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"


class PersonIndexGPU : public PersonIndex {
DISALLOW_COPY_AND_ASSIGN(PersonIndexGPU)

public:
    PROPERTY_REF(ThrustTVectorHost<Person*>, h_persons);
    PROPERTY_REF(ThrustTVectorHost<glm::mat4>, h_person_models);
    PROPERTY_REF(ThrustTVectorHost<glm::vec4>, h_person_colors);
    PROPERTY_REF(ThrustTVectorHost<long>, h_person_ids);
    PROPERTY_REF(ThrustTVectorHost<int>, h_person_host_states);
    PROPERTY_REF(ThrustTVectorHost<int>, h_person_ages);
    PROPERTY_REF(ThrustTVectorHost<int>, h_person_age_classes);
    PROPERTY_REF(ThrustTVectorHost<int>, h_person_locations);
    PROPERTY_REF(ThrustTVectorHost<int>, h_person_residence_locations);

    PROPERTY_REF(ThrustTVectorHost<double>, h_person_innate_relative_biting_rates);
    PROPERTY_REF(ThrustTVectorHost<double>, h_person_current_relative_biting_rates);

    //on DEVICE CUDA
    PROPERTY_REF(ThrustTVectorDevice<glm::mat4>,buffer_person_models);
    PROPERTY_REF(ThrustTVectorDevice<glm::vec4>,buffer_person_colors);

public:
    PersonIndexGPU();

    virtual ~PersonIndexGPU();

void add(Person *p);

void add(Person *p, const int &location, const Person::HostStates &host_state, const int &age_class);

virtual void remove(Person *p);

virtual std::size_t size() const;

virtual void update();

void remove_without_set_index(Person *p);

virtual void notify_change(Person *p, const Person::Property &property, const void *oldValue, const void *newValue);

void change_property(Person *p, const int &location, const Person::HostStates &host_state, const int &age_class);

private:
};


#endif //MASS_PERSONINDEXGPU_H
