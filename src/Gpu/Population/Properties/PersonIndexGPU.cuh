//
// Created by kient on 12/9/2023.
//

#ifndef PERSONINDEXGPU_CUH
#define PERSONINDEXGPU_CUH

#include "PersonIndex.cuh"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

namespace GPU{
    class PersonIndexGPU;
    class Person;
}

class GPU::PersonIndexGPU : public GPU::PersonIndex {
DISALLOW_COPY_AND_ASSIGN(PersonIndexGPU)

public:
    PROPERTY_REF(TVector<GPU::Person*>, h_persons);
    PROPERTY_REF(TVector<glm::mat4>, h_person_models);
    PROPERTY_REF(TVector<glm::vec4>, h_person_colors);
    PROPERTY_REF(TVector<long>, h_person_ids);
    PROPERTY_REF(TVector<int>, h_person_host_states);
    PROPERTY_REF(TVector<int>, h_person_ages);
    PROPERTY_REF(TVector<int>, h_person_age_classes);
    PROPERTY_REF(TVector<int>, h_person_locations);
    PROPERTY_REF(TVector<int>, h_person_residence_locations);
    PROPERTY_REF(TVector<int>, h_person_n_parasites);
    PROPERTY_REF(TVector<double>, h_person_log10_total_infectious_denstiy);

    PROPERTY_REF(TVector<double>, h_person_innate_relative_biting_rates);
    PROPERTY_REF(TVector<double>, h_person_current_relative_biting_rates);
    PROPERTY_REF(TVector<GPU::SingleHostClonalParasitePopulations*>, h_person_parasites);


    //on DEVICE CUDA
    PROPERTY_REF(ThrustTVectorDevice<glm::mat4>,buffer_person_models);
    PROPERTY_REF(ThrustTVectorDevice<glm::vec4>,buffer_person_colors);

public:
    PersonIndexGPU();

    virtual ~PersonIndexGPU();

void add(GPU::Person *p);

void add(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);

virtual void remove(GPU::Person *p);

virtual std::size_t size() const;

virtual void update();

void remove_without_set_index(GPU::Person *p);

virtual void notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue);

void change_property(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);

private:
};


#endif //PERSONINDEXGPU_H
