//
// Created by kient on 12/9/2023.
//

#ifndef PERSONINDEXGPU_CUH
#define PERSONINDEXGPU_CUH

#include "PersonIndex.cuh"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

#define MAX_PARASITE_PER_PERSON 10
#define MAX_DRUG_NUMBER 8
#define MAX_DRUG_PER_ACT 3
#define MAX_GENOTYPE_LOCUS 38

namespace GPU{
    class PersonIndexGPU;
    class Person;

    /*
     * Using MAX_PARASITE_PER_PERSON and MAX_DRUG_PER_ACT to define the size of the array
     * This is to remove overhead of using vector
     * Update the size of the array if needed
     * */
    struct PersonUpdateInfo{
        long person_id;
        int person_index;
        int person_latest_update_time;
        double person_latest_immune_value;
        /* for parasite update */
        int parasite_size;
        long parasite_id;
        int parasite_current_index;
        int parasite_update_function_type[MAX_PARASITE_PER_PERSON];
        double parasite_last_update_log10_parasite_density[MAX_PARASITE_PER_PERSON];
        double parasite_genotype_fitness_multiple_infection[MAX_PARASITE_PER_PERSON];
        double parasite_gametocyte_level[MAX_PARASITE_PER_PERSON];
        /* for drug update */
        int drug_in_blood_size;
        double drug_starting_value[MAX_DRUG_NUMBER];
        int drug_dosing_days[MAX_DRUG_NUMBER];
        double drug_last_update_value[MAX_DRUG_NUMBER];
        int drug_last_update_time[MAX_DRUG_NUMBER];
        int drug_start_time[MAX_DRUG_NUMBER];
        int drug_end_time[MAX_DRUG_NUMBER];
        double drug_half_life[MAX_DRUG_NUMBER];
        double drug_rand_uniform_1[MAX_DRUG_NUMBER];
        double drug_rand_uniform_2[MAX_DRUG_NUMBER];
        int drug_in_blood_type_id[MAX_DRUG_PER_ACT];
        double drug_epsilon;
        /* for parasite update with drug */
        char parasite_genotype[MAX_PARASITE_PER_PERSON][MAX_GENOTYPE_LOCUS+1];
    };
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
    PROPERTY_REF(TVector<GPU::PersonUpdateInfo>, h_person_update_info);

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
    GPU::PersonUpdateInfo init_person_update_info(GPU::Person *p);
};


#endif //PERSONINDEXGPU_H
