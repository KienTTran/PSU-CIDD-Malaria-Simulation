/*
 * File:   TypeDef.h
 * Author: nguyentran
 *
 * Created on April 17, 2013, 10:17 AM
 */

#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <array>
#include <list>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include <thrust/tuple.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>
#include <GL/glew.h>

namespace GPU{
    class Person;
    class PersonIndex;
    class IStrategy;
    class Event;
    class Therapy;
    class Drug;
    class PersonUpdateInfo;
};

class Person;

class PersonIndex;

class Event;

class Reporter;

class Drug;

class IStrategy;

class Therapy;


template <typename T>
using ThrustTuple = thrust::tuple<T>;
template <typename T>
using ThrustTupleVectorDevice = thrust::device_vector<ThrustTuple<T>>;
template <typename T>
using ThrustTupleVector = std::vector<ThrustTuple<T>>;

template <typename T, typename T2>
using ThrustTuple2 = thrust::tuple<T, T2>;
template <typename T, typename T2>
using ThrustTuple2Vector = std::vector<ThrustTuple2<T,T2>>;
template <typename T, typename T2>
using ThrustTuple2VectorDevice = thrust::device_vector<ThrustTuple2<T,T2>>;

template <typename T, typename T2, typename T3>
using ThrustTuple3 = thrust::tuple<T, T2, T3>;
template <typename T, typename T2, typename T3>
using ThrustTuple3Vector=  std::vector<ThrustTuple3<T,T2,T3>>;
template <typename T, typename T2, typename T3>
using ThrustTuple3VectorDevice = thrust::device_vector<ThrustTuple3<T,T2,T3>>;

template <typename T, typename T2, typename T3, typename T4>
using ThrustTuple4 = thrust::tuple<T, T2, T3, T4>;
template <typename T, typename T2, typename T3, typename T4>
using ThrustTuple4Vector = std::vector<ThrustTuple4<T,T2,T3,T4>>;
template <typename T, typename T2, typename T3, typename T4>
using ThrustTuple4VectorDevice = thrust::device_vector<ThrustTuple4<T,T2,T3,T4>>;

template <typename T, typename T2, typename T3, typename T4, typename T5>
using ThrustTuple5 = thrust::tuple<T, T2, T3, T4, T5>;
template <typename T, typename T2, typename T3, typename T4, typename T5>
using ThrustTuple5VectorHost = thrust::host_vector<ThrustTuple5<T,T2,T3,T4,T5>>;
template <typename T, typename T2, typename T3, typename T4, typename T5>
using ThrustTuple5VectorDevice = thrust::device_vector<ThrustTuple5<T,T2,T3,T4,T5>>;

template <typename T>
using TVector = std::vector<T>;
template <typename T>
using TVector2 = std::vector<TVector<T>>;
template <typename T>
using TVector3 = std::vector<TVector2<T>>;
template <typename T>
using TVector4 = std::vector<TVector3<T>>;
template <typename T>
using ThrustTVectorDevice = thrust::device_vector<T>;

typedef unsigned long ul;

typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleVector2;
typedef std::vector<DoubleVector2> DoubleVector3;
typedef std::vector<int> IntVector;
typedef std::vector<int> *IntVectorPtr;
typedef std::vector<IntVector> IntVector2;
typedef std::vector<IntVector2> IntVector3;
typedef std::vector<IntVector *> IntVectorPtrVector;
typedef std::vector<IntVector> *IntVector2Ptr;
typedef std::vector<unsigned int> UIntVector;

typedef std::vector<ul> LongVector;
typedef std::vector<LongVector> LongVector2;

typedef std::vector<std::string> StringVector;
typedef std::vector<StringVector> StringVector2;

typedef std::map<int, int> IntIntMap;

typedef std::vector<Person *> PersonPtrVector;
typedef std::vector<PersonPtrVector> PersonPtrVector2;
typedef std::vector<PersonPtrVector2> PersonPtrVector3;
typedef std::vector<PersonPtrVector3> PersonPtrVector4;

typedef std::vector<GPU::Person *> GPUPersonPtrVector;
typedef std::vector<GPUPersonPtrVector> GPUPersonPtrVector2;
typedef std::vector<GPUPersonPtrVector2> GPUPersonPtrVector3;
typedef std::vector<GPUPersonPtrVector3> GPUPersonPtrVector4;

typedef std::vector<Event *> EventPtrVector;
typedef std::vector<EventPtrVector> EventPtrVector2;

typedef std::vector<GPU::Event *> GPUEventPtrVector;
typedef std::vector<GPUEventPtrVector> GPUEventPtrVector2;

typedef std::vector<Reporter *> ReporterPtrVector;

typedef std::list<PersonIndex *> PersonIndexPtrList;
typedef std::list<GPU::PersonIndex *> GPUPersonIndexPtrList;

typedef std::map<int, Drug *> DrugPtrMap;

typedef std::vector<Therapy *> TherapyPtrVector;
typedef std::vector<IStrategy *> StrategyPtrVector;

typedef std::map<int, GPU::Drug *> GPUDrugPtrMap;
typedef std::vector<GPU::Therapy *> GPUTherapyPtrVector;
typedef std::vector<GPU::IStrategy *> GPUStrategyPtrVector;

template <typename T>
using TVector = std::vector<T>;
template <typename T>
using TVector2 = std::vector<TVector<T>>;
template <typename T>
using TVector3 = std::vector<TVector2<T>>;

using GenotypeDBInfo = thrust::tuple<int,char*,double>;
template<class T>
using ThrustTVectorHostPinned = thrust::host_vector<T, thrust::mr::stateless_resource_allocator<T,
        thrust::system::cuda::universal_host_pinned_memory_resource>>;

struct GPUConfig{
    int n_threads {1024};
    int n_streams {8};
    int n_people_1_batch {1000000};
    double pre_allocated_mem_ratio {1.0};
    friend std::ostream &operator<<(std::ostream &os, const GPUConfig &mcf) {
        os << "gpu_config";
        return os;
    }
};

struct RenderConfig{
    int window_width;
    int window_height;
    bool display_gui;
    bool display_plot;
    bool close_window_on_finish;
    double point_coord;
    std::string vertex_shader_path;
    std::string fragment_shader_path;
    friend std::ostream &operator<<(std::ostream &os, const RenderConfig &mcf) {
        os << "render_config";
        return os;
    }
};

struct DebugConfig{
    int width;
    int height;
    bool enable_random_position_update;
    bool enable_debug_text;
    bool enable_debug_render;
    bool enable_debug_render_text;
    int log_interval;
    friend std::ostream &operator<<(std::ostream &os, const DebugConfig &mcf) {
        os << "debug_config";
        return os;
    }
};

struct RasterDb {
    friend std::ostream &operator<<(std::ostream &os, const RasterDb &rdb) {
        os << "raster_db";
        return os;
    }
};

struct SeasonalInfo {
  bool enable { false };
  DoubleVector A;
  DoubleVector B;
  DoubleVector C;
  DoubleVector phi;
  DoubleVector min_value;

  friend std::ostream &operator<<(std::ostream &os, const SeasonalInfo &seasonal_info);
};

inline std::ostream &operator<<(std::ostream &os, const SeasonalInfo &seasonal_info) {
  os << "seasonal_info: ";
  return os;
}

struct ImmuneSystemInformation {
  double acquire_rate { -1 };
  /*
   * Use array instead of vector to use in GPU
   * Check void immune_system_information::set_value(const YAML::Node &node)
   * when size is different from 81
   *
   * */
  double acquire_rate_by_age[81];
//  std::vector<double> acquire_rate_by_age;
  double decay_rate { -1 };

  double duration_for_fully_immune { -1 };
  double duration_for_naive { -1 };

  //    double mean_initial_condition;
  //    double sd_initial_condition;

  double immune_inflation_rate { -1 };

  double min_clinical_probability { -1 };
  double max_clinical_probability { -1 };

  double immune_effect_on_progression_to_clinical { -1 };

  double c_min { -1 };
  double c_max { -1 };

  double alpha_immune { -1 };
  double beta_immune { -1 };

  double age_mature_immunity { -1 };
  double factor_effect_age_mature_immunity { -1 };
};

/*
 * If there is a vector or array in this struct,
 * please allocate it to GPU device before using in population.cu
 * */
struct ParasiteDensityLevel {
  double log_parasite_density_cured;
  double log_parasite_density_from_liver;
  double log_parasite_density_asymptomatic;
  double log_parasite_density_clinical;
  double log_parasite_density_clinical_from;
  double log_parasite_density_clinical_to;
  double log_parasite_density_detectable;
  double log_parasite_density_detectable_pfpr;
  double log_parasite_density_pyrogenic;

  friend std::ostream &operator<<(std::ostream &os, const ParasiteDensityLevel &pdl) {
    os << "[" << pdl.log_parasite_density_cured << "," << pdl.log_parasite_density_from_liver << ","
       << pdl.log_parasite_density_asymptomatic << "," << pdl.log_parasite_density_clinical << ","
       << pdl.log_parasite_density_clinical_from << "," << pdl.log_parasite_density_clinical_to << ","
       << pdl.log_parasite_density_detectable << "," << pdl.log_parasite_density_detectable_pfpr << ","
       << pdl.log_parasite_density_pyrogenic << "]";
    return os;
  }
};

struct RelativeBittingInformation {
  double max_relative_biting_value;
  double min_relative_biting_value;
  int number_of_biting_levels;

  double scale;

  double mean;
  double sd;

  double gamma_a;
  double gamma_b;
};

struct RelativeMovingInformation {
  double max_relative_moving_value;
  int number_of_moving_levels;

  //  biting_level_distribution:
  //  #  distribution: Exponential
  //    distribution: Gamma
  //    Exponential:
  double scale;

  double mean;
  double sd;
  DoubleVector v_moving_level_value;
  DoubleVector v_moving_level_density;

  double circulation_percent;
  double length_of_stay_mean;
  double length_of_stay_sd;
  double length_of_stay_theta;
  double length_of_stay_k;
};

struct InitialParasiteInfo {
  int location;
  int parasite_type_id;
  double prevalence;

  InitialParasiteInfo() : location(-1), parasite_type_id(-1), prevalence(-1.0) {};

  InitialParasiteInfo(const int loc, const int p_type, const double pre)
      : location(loc), parasite_type_id(p_type), prevalence(pre) {};
};

struct RelativeInfectivity {
  double sigma;
  double ro_star;

  friend std::ostream &operator<<(std::ostream &os, const RelativeInfectivity &e) {
    os << "[" << e.sigma << "," << e.ro_star << "]";
    return os;
  }
};

struct AaPositionInfo {
  int position { -1 };
  std::vector<char> amino_acids;
  std::vector<double> daily_crs;
  std::map<int, std::vector<double>> multiplicative_effect_on_EC50;
  friend std::ostream &operator<<(std::ostream &os, const AaPositionInfo &aa) { return os; }
};

struct GeneInfo {
  std::string name;
  int chromosome { -1 };
  int max_copies { 1 };
  std::vector<double> cnv_daily_crs;
  std::map<int, std::vector<double>> cnv_multiplicative_effect_on_EC50;
  std::vector<AaPositionInfo> aa_position_infos;
  std::map<int, double> multiplicative_effect_on_EC50_for_2_or_more_mutations;
  friend std::ostream &operator<<(std::ostream &os, const GeneInfo &aa) { return os; }
  double average_daily_crs { -1.0 };
};

struct ChromosomeInfo {
  std::vector<GeneInfo> gene_infos;

  friend std::ostream &operator<<(std::ostream &os, const ChromosomeInfo &chromosome) { return os; }
};
struct PfGeneInfo {
  std::array<ChromosomeInfo, 14> chromosome_infos {};

  friend std::ostream &operator<<(std::ostream &os, const PfGeneInfo &geneInfo) { return os; }
  int calculate_aa_pos(int chromosome_id, int gene_id, int aa_pos_id, int aa_id) {
    auto result = 0;
    for (int ii = 0; ii < chromosome_id; ++ii) {
      result += chromosome_infos[ii].gene_infos.size() > 1 ? chromosome_infos[ii].gene_infos.size() - 1 : 0;  // for ','
      for (auto &gene_info : chromosome_infos[ii].gene_infos) {
        result += gene_info.aa_position_infos.size();
        result += gene_info.max_copies > 1 ? 1 : 0;  // for copy number
      }
    }
    result += chromosome_id;  // for "|"

    // final chromosome
    for (int ii = 0; ii < gene_id; ++ii) {
      result += chromosome_infos[chromosome_id].gene_infos[ii].aa_position_infos.size();
      result += chromosome_infos[chromosome_id].gene_infos[ii].max_copies > 1 ? 1 : 0;  // for copy number
      result += 1;                                                                      // for ","
    }

    // final gene
    result += aa_id;
    return result;
  }
};

struct OverrideEC50Pattern {
  std::string pattern;
  int drug_id;
  double ec50;

  friend std::ostream &operator<<(std::ostream &os, const OverrideEC50Pattern &pattern) {
    os << "[";
    os << pattern.pattern << ", ";
    os << pattern.drug_id << ", ";
    os << pattern.ec50;
    return os << "]";
  }
};

typedef std::vector<OverrideEC50Pattern> OverrideEC50Patterns;

struct MosquitoConfig{
  std::vector<double> interrupted_feeding_rate;
  std::string interrupted_feeding_rate_raster {""};
  int prmc_size {20};
  friend std::ostream &operator<<(std::ostream &os, const MosquitoConfig &mcf) {
    for(auto rate : mcf.interrupted_feeding_rate){
      os << "Mosquito size: " << mcf.prmc_size << ", IF rate: " << rate;
    }
    return os;
  }
};

#endif /* TYPEDEF_H */
