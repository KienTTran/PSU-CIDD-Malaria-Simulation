//
// Created by kient on 12/31/2023.
//

#ifndef GPUPOPULATIONKERNEL_CUH
#define GPUPOPULATIONKERNEL_CUH

#include "Core/TypeDef.h"
#include "Gpu/Therapies/DrugType.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"

namespace GPU{
    class PopulationKernel;
}

class GPU::PopulationKernel {
public:
    /*
     * Variables for circulations
     * */
    ThrustTVectorDevice<int> d_ce_popsize_residence_by_location;
    ThrustTVectorDevice<double> d_ce_spatial_model_parameters;
    ThrustTVectorDevice<double> d_ce_spatial_model_travels;
    ThrustTVectorDevice<int> d_ce_spatial_districts;
    ThrustTVectorDevice<double> d_ce_spatial_distances;
    ThrustTVectorDevice<int> d_ce_popsize_by_moving_level;
    ThrustTVectorDevice<double> d_ce_moving_level_value;
    ThrustTVectorDevice<int> d_ce_n_circulations_from_locations;
    ThrustTVectorDevice<double> d_ce_relative_outmovement_from_target;
    ThrustTVectorDevice<int> d_ce_all_location_from_indices;
    ThrustTVectorDevice<int> d_ce_all_location_target_indices;
    ThrustTVectorDevice<int> d_ce_all_moving_levels;

    /*
     * Variables for update individuals
     * */
    TVector<GPU::DrugType::ResistantAALocation> h_drug_res_aa_loc;
    ThrustTVectorDevice<GPU::DrugType::ResistantAALocation> d_drug_res_aa_loc;
    TVector<int> h_drug_res_aa_loc_index;
    ThrustTVectorDevice<int> d_drug_res_aa_loc_index;
    TVector<int> h_gen_gene_size;
    ThrustTVectorDevice<int> d_gen_gene_size;
    TVector<int> h_gen_max_copies;
    ThrustTVectorDevice<int> d_gen_max_copies;
    TVector<int> h_gen_aa_size;
    ThrustTVectorDevice<int> d_gen_aa_size;
    TVector<std::string> h_gen_aa_test;
    TVector<int> h_gen_aa_int;
    ThrustTVectorDevice<int> d_gen_aa_int;
    TVector<int> h_gen_aa_int_start_index;
    ThrustTVectorDevice<int> d_gen_aa_int_start_index;
    char* d_gen_mutation_mask{};
    ImmuneSystemInformation *d_immune_system_information{};
    GPU::PersonIndexGPU *pi{};
    GPU::PersonUpdateInfo *d_buffer_person_update_info_stream{};
    GPU::Person *d_buffer_person_stream{};
    cudaStream_t *d_streams{};
    cudaEvent_t start_event{}, stop_event{};
    float elapsed_time{};

    /*
     * Variables for update FOI
     * */
    /* This vector store relative biting, moving and foi of all people by location */
    TVector<ThrustTuple4<int,double,double,double>> h_sum_biting_moving_foi_by_loc;
    ThrustTVectorDevice<ThrustTuple4<int,double,double,double>> d_sum_biting_moving_foi_by_loc;
    TVector<TVector<double>> individual_foi_by_location;
    TVector<TVector<double>> individual_relative_biting_by_location;
    TVector<TVector<double>> individual_relative_moving_by_location;
    TVector<TVector<double>> force_of_infection_for_N_days_by_location;

public:
    PopulationKernel();
    ~PopulationKernel() = default;
public:
    void init();

    /*
     * for circulations
     * */
    void calculate_circulate_locations(int n_locations,
                                       ThrustTVectorDevice<int> &d_n_circulations_by_location,
                                       ThrustTVectorDevice<double> &d_relative_outmovement_from_to,
                                       ThrustTVectorDevice<int> &d_all_location_from_indices,
                                       ThrustTVectorDevice<int> &d_all_location_to_indices);
    void calculate_moving_level_density(ThrustTuple2VectorDevice<int,int> d_circulation_indices,
                                        ThrustTVectorDevice<double> &d_moving_level_density);
    void perform_circulation_event();

    /*
     * for infection
     * */
//    void calculate_n_person_bitten_today(int n_locations,
//                                          ThrustTVectorDevice<double> &d_foi_all_locations,
//                                          ThrustTVectorDevice<int> &d_n_person_bitten_today_all_locations);
//    void perform_infection_event();
    void update_all_individuals();
    void update_current_foi();
    void persist_current_force_of_infection_to_use_N_days_later();
};



#endif //GPUPOPULATION_CUH
