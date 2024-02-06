//
// Created by kient on 12/9/2023.
//

#ifndef PERSONINDEXGPU_CUH
#define PERSONINDEXGPU_CUH

#include "PersonIndex.cuh"
#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"

#define MAX_PARASITE_PER_PERSON 10
#define MAX_DRUG_NUMBER 9
#define MAX_DRUG_PER_ACT 3
#define MAX_GENOTYPE_LOCUS 38



namespace GPU{
    class PersonIndexGPU;
    class Person;

    /*
     * Using MAX_PARASITE_PER_PERSON and MAX_DRUG_PER_ACT to define the size of the array
     * This is to remove overhead of using vector in cuda kernel
     * Update the size of the array if needed
     * */
    struct PersonUpdateInfo{
        long person_id;
        int person_index;
        int person_age;
        int person_age_class;
        int person_host_state;
        int person_location;
        int person_residence_location;
        int person_birthday;
        int person_moving_level;
        double person_innate_relative_biting_rate;
        double person_current_relative_biting_rate;
        bool person_using_age_dependent_biting_level;
        int person_latest_update_time;
        char person_liver_parasite_genotype[MAX_GENOTYPE_LOCUS+1];
        bool person_has_drug_in_blood;
        /*
         * for parasite update
         * parasite_ is clonal parasite variable
         * parasites_ is single host clonal parasite populations variable
         * one single host clonal parasite populations can have MAX_PARASITE_PER_PERSON clonal parasites
         * */
        double limit_epsilon;
        double LOG_ZERO_PARASITE_DENSITY;
        int parasites_size;
        int parasites_current_index;
        double parasites_log10_total_infectious_density;
        int parasites_genotype_mutated_number;
        long parasite_id[MAX_PARASITE_PER_PERSON];
        int parasite_update_function_type[MAX_PARASITE_PER_PERSON];
        double parasite_last_update_log10_parasite_density[MAX_PARASITE_PER_PERSON];
        double parasite_genotype_fitness_multiple_infection[MAX_PARASITE_PER_PERSON];
        double parasite_gametocyte_level[MAX_PARASITE_PER_PERSON];
        double parasite_log10_infectious_density[MAX_PARASITE_PER_PERSON];
        /* for parasite update with drug */
        char parasite_genotype[MAX_PARASITE_PER_PERSON][MAX_GENOTYPE_LOCUS+1];
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
        /* for immune system update */
        double immune_system_component_latest_value;
        int immune_system_component_type;
        bool immune_system_is_increased;
    };
}

class GPU::PersonIndexGPU : public GPU::PersonIndex {
DISALLOW_COPY_AND_ASSIGN(PersonIndexGPU)

public:
    PROPERTY_REF(TVector<GPU::Person*>, h_persons);
    PROPERTY_REF(TVector<glm::mat4>, h_person_models);
    PROPERTY_REF(TVector<glm::vec4>, h_person_colors);
    PROPERTY_REF(TVector<int>, h_person_residence_locations);
    PROPERTY_REF(TVector<GPU::PersonUpdateInfo>, h_person_update_info);

    //on DEVICE CUDA
    PROPERTY_REF(ThrustTVectorDevice<glm::mat4>,buffer_person_models);
    PROPERTY_REF(ThrustTVectorDevice<glm::vec4>,buffer_person_colors);

public:
    PersonIndexGPU();

    virtual ~PersonIndexGPU();

void add(GPU::Person *p);

//void add(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);

virtual void remove(GPU::Person *p);

virtual std::size_t size() const;

virtual void update();

void remove_without_set_index(GPU::Person *p);

virtual void notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue, const void *newValue);

//void change_property(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state, const int &age_class);

public:
    void update_person();

private:
    GPU::PersonUpdateInfo init_person_update_info(GPU::Person *p);

};


#endif //PERSONINDEXGPU_H
