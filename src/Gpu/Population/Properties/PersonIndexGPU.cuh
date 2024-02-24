//
// Created by kient on 12/9/2023.
//

#ifndef PERSONINDEXGPU_CUH
#define PERSONINDEXGPU_CUH

#include "PersonIndex.cuh"
#include "Gpu/Therapies/Drug.cuh"
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
     * Using MAX_PARASITE_PER_PERSON, MAX_DRUG_NUMBERand MAX_DRUG_PER_ACT to define the size of
     * drug and parasite array in the struct PersonUpdateInfo
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
        double person_current_relative_moving_rate;
        double person_current_relative_infectivity;
        double person_current_foi;
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
        char parasite_genotype_current[MAX_PARASITE_PER_PERSON][MAX_GENOTYPE_LOCUS+1];
        char parasite_genotype_new[MAX_PARASITE_PER_PERSON][MAX_GENOTYPE_LOCUS+1];
        double parasite_genotype_daily_fitness_multiple_infection[MAX_PARASITE_PER_PERSON];
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
        int drug_in_blood_type_id_current_index;
        /* for immune system update */
        double immune_system_component_latest_value;
        int immune_system_component_type;
        bool immune_system_is_increased;

        void add_parasite_to_blood(const int &id_,
                                   const int &index_,
                                   const int &latest_update_time_,
                                   GPU::SingleHostClonalParasitePopulations* all_clonal_parasite_populations_,
                                   GPU::ClonalParasitePopulation* blood_parasite){
          person_id = id_;
          person_index = index_;
          person_latest_update_time = latest_update_time_;
          parasites_size = all_clonal_parasite_populations_->size();
          parasites_current_index = blood_parasite->index();
          parasite_id[blood_parasite->index()] = blood_parasite->id();
        }

        void add_parasite_genotype(GPU::Genotype* parasite_type, int parasite_index){
          /* Convert string to char[] with NULL terminated */
          std::copy(parasite_type->aa_sequence.begin(),
                    parasite_type->aa_sequence.end(),
                    parasite_genotype_current[parasite_index]);
          parasite_genotype_current[parasite_index][MAX_GENOTYPE_LOCUS] = '\0';
          std::copy(parasite_type->aa_sequence.begin(),
                    parasite_type->aa_sequence.end(),
                    parasite_genotype_new[parasite_index]);
          parasite_genotype_new[parasite_index][MAX_GENOTYPE_LOCUS] = '\0';
          parasite_genotype_fitness_multiple_infection[parasite_index] = parasite_type->daily_fitness_multiple_infection;
        }

        void add_drug_in_blood(GPUDrugPtrMap *drugs, const int &drug_id){
            if(drug_in_blood_type_id_current_index >= MAX_DRUG_PER_ACT || drug_in_blood_size >= MAX_DRUG_NUMBER){
                LOG(FATAL) << fmt::format("Exceed maximum drug per act ({})", MAX_DRUG_PER_ACT);
            }
            drug_dosing_days[drug_id] = drugs->at(drug_id)->dosing_days();
            drug_half_life[drug_id] = drugs->at(drug_id)->drug_type()->drug_half_life();
            drug_starting_value[drug_id] = drugs->at(drug_id)->starting_value();
            drug_last_update_value[drug_id] = drugs->at(drug_id)->last_update_value();
            drug_start_time[drug_id] = drugs->at(drug_id)->start_time();
            drug_end_time[drug_id] = drugs->at(drug_id)->end_time();
            drug_rand_uniform_1[drug_id] = 0.1;//Model::RANDOM->random_uniform_double(-0.2, 0.2),
            drug_rand_uniform_2[drug_id] = 0.1;//Model::RANDOM->random_uniform_double(0, 0.1))
            drug_last_update_time[drug_id] = drugs->at(drug_id)->last_update_time();
            drug_in_blood_type_id[drug_in_blood_type_id_current_index] = drug_id;
            drug_in_blood_type_id_current_index++;
            drug_in_blood_size++;
            person_has_drug_in_blood = true;
        }
        void remove_drug_in_blood(const int &index){
            if(drug_in_blood_type_id_current_index <= 0 || drug_in_blood_size <= 0){
                LOG(FATAL) << "Drug in blood is empty";
            }
            for(int i = index; i < drug_in_blood_type_id_current_index - 1; i++){
                drug_dosing_days[i] = drug_dosing_days[i+1];
                drug_half_life[i] = drug_half_life[i+1];
                drug_starting_value[i] = drug_starting_value[i+1];
                drug_last_update_value[i] = drug_last_update_value[i+1];
                drug_start_time[i] = drug_start_time[i+1];
                drug_end_time[i] = drug_end_time[i+1];
                drug_rand_uniform_1[i] = drug_rand_uniform_1[i+1];
                drug_rand_uniform_2[i] = drug_rand_uniform_2[i+1];
                drug_last_update_time[i] = drug_last_update_time[i+1];
                drug_in_blood_type_id[i] = drug_in_blood_type_id[i+1];
            }
            drug_in_blood_type_id_current_index--;
            drug_in_blood_size--;
            if(drug_in_blood_type_id_current_index == 0 && drug_in_blood_size == 0){
                person_has_drug_in_blood = false;
            }
        }
    };
}


class GPU::PersonIndexGPU : public GPU::PersonIndex {
DISALLOW_COPY_AND_ASSIGN(PersonIndexGPU)

public:
    PROPERTY_REF(TVector<GPU::Person*>, h_persons);
    PROPERTY_REF(ThrustTVectorHostPinned<glm::mat4>, h_person_models);
    PROPERTY_REF(ThrustTVectorHostPinned<glm::vec4>, h_person_colors);
    PROPERTY_REF(TVector<int>, h_person_residence_locations);
    PROPERTY_REF(ThrustTVectorHostPinned<GPU::PersonUpdateInfo>, h_person_update_info);
    PROPERTY_REF(TVector<GenotypeDBInfo>, h_genotype_db);
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
