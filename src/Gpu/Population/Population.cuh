//
// Created by kient on 12/31/2023.
//

#ifndef MASS_GPUPOPULATION_CUH
#define MASS_GPUPOPULATION_CUH


#include "Core/TypeDef.h"

namespace GPU{
    class Population;
}

class GPU::Population {
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
     * Variables for infection
     * */
    TVector<double> h_ie_foi_N_days_all_locations;/* index is loc_index*n_tracking_day+track_day_index */
    ThrustTVectorDevice<int> d_ie_foi_N_days_all_locations;
public:
    Population();
    ~Population() = default;
public:
    void init();

    void update_current_foi();
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
    void calculate_n_person_bitten_today(int n_locations,
                                          ThrustTVectorDevice<double> &d_foi_all_locations,
                                          ThrustTVectorDevice<int> &d_n_person_bitten_today_all_locations);
    void perform_infection_event();
};


#endif //MASS_GPUPOPULATION_CUH
