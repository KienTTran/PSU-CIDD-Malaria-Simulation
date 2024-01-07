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
private:
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
    ThrustTVectorDevice<int> d_infection_popsize_by_location;
public:
    Population();
    ~Population() = default;
public:
    void init();
    void calculate_circulate_locations(int n_locations,ThrustTVectorDevice<int> &d_n_circulations_by_location,ThrustTVectorDevice<double> &d_relative_outmovement_from_to,
                            ThrustTVectorDevice<int> &d_all_location_from_indices,ThrustTVectorDevice<int> &d_all_location_to_indices);
    void calculate_moving_level_density(ThrustT2TupleVectorDevice<int,int> d_circulation_indices,ThrustTVectorDevice<double> &d_moving_level_density);
    void perform_circulation_event();
    void perform_infection_event();
};


#endif //MASS_GPUPOPULATION_CUH
