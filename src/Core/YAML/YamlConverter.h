#ifndef YAMLCONVERTER_H
#define YAMLCONVERTER_H

#include <date/date.h>
#include <yaml-cpp/yaml.h>

#include <cmath>

#include "Core/TypeDef.h"
#include "GIS/SpatialData.h"
#include "Spatial/Location.h"
#include "easylogging++.h"

namespace YAML {
template <>
struct convert<date::sys_days> {
  static Node encode(const date::sys_days& rhs) {
    Node node;
    node.push_back(date::format("%Y/%m/%d", rhs));
    return node;
  }

  static bool decode(const Node& node, date::sys_days& rhs) {
    if (!node.IsScalar()) {
      return false;
    }
    std::stringstream ss(node.as<std::string>());
    date::from_stream(ss, "%Y/%m/%d", rhs);
    return true;
  }
};

template <>
struct convert<date::year_month_day> {
  static Node encode(const date::year_month_day& rhs) {
    Node node;
    node.push_back(date::format("%Y/%m/%d", rhs));
    return node;
  }

  static bool decode(const Node& node, date::year_month_day& rhs) {
    if (!node.IsScalar()) {
      return false;
    }
    std::stringstream ss(node.as<std::string>());
    from_stream(ss, "%Y/%m/%d", rhs);
    return true;
  }
};

template <>
struct convert<GPUConfig> {
    static Node encode(const GPUConfig& gcfe) {
        Node node;
        node.push_back("gpu_config");
        return node;
    }
    static bool decode(const Node& node, GPUConfig& gcfd) {
        gcfd.n_threads = node["n_threads"].as<int>();
        gcfd.people_1_batch = node["people_1_batch"].as<int>();
        gcfd.pre_allocated_mem_ratio = node["pre_allocated_mem_ratio"].as<double>();
        gcfd.n_parasites_per_person = node["n_parasites_per_person"].as<int>();
        return true;
    }
};


template <>
struct convert<RenderConfig> {
    static Node encode(const RenderConfig& rcfe) {
        Node node;
        node.push_back("render_config");
        return node;
    }
    static bool decode(const Node& node, RenderConfig& rcfd) {
        rcfd.display_gui = node["display_gui"].as<bool>();
        rcfd.close_window_on_finish = node["close_window_on_finish"].as<bool>();
        rcfd.point_coord = node["point_coord"].as<double>();
        rcfd.window_width = node["window_width"].as<int>();
        rcfd.window_height = node["window_height"].as<int>();
        rcfd.fragment_shader_path = node["shaders"]["fragment"].as<std::string>();
        rcfd.vertex_shader_path = node["shaders"]["vertex"].as<std::string>();
        return true;
    }
};

template <>
struct convert<DebugConfig> {
    static Node encode(const DebugConfig& dgbcfe) {
        Node node;
        node.push_back("debug_config");
        return node;
    }
    static bool decode(const Node& node, DebugConfig& dgbcfd) {
        dgbcfd.enable_debug_render = node["enable_debug_render"].as<bool>();
        dgbcfd.enable_debug_text = node["enable_debug_text"].as<bool>();
        dgbcfd.enable_debug_render_text = node["enable_debug_render_text"].as<bool>();
        dgbcfd.enable_update = node["enable_update"].as<bool>();
        dgbcfd.height = node["height"].as<int>();
        dgbcfd.width = node["width"].as<int>();
        dgbcfd.width = node["width"].as<int>();
        dgbcfd.log_interval = node["log_interval"].as<int>();
        return true;
    }
};

template<>
struct convert<RasterDb> {
    static Node encode(const RasterDb &rdb) {
      Node node;
      node.push_back("raster_db");
      return node;
    }
    static bool decode(const Node &node, RasterDb &rdb) {
      return SpatialData::get_instance().parse(node);
    }
};

template<>
struct convert<std::vector<Spatial::Location>> {
    static Node encode(const std::vector<Spatial::Location> &rhs) {
        Node node;
        node.push_back("location_db");
        return node;
    }

    // Decode the contents of the location_db node
    static bool decode(const Node &node, std::vector<Spatial::Location> &location_db) {

        // Check to see if the location informatoin has already been entered, this implies
        // that there was a raster_db node present and an error starte exists
        if (location_db.size() != 0) {
            throw std::runtime_error("location_db has already been instantiated, is a raster_db present in the file?");
        }

        // If the user is supplying raster data, location_info will likely not be there.
        // Since we are just decoding the file, just focus on loading the data and defer
        // validation of the file to the recipent of the data
        auto number_of_locations = node["location_info"].size();
        for (std::size_t i = 0; i < number_of_locations; i++) {
            location_db.emplace_back(node["location_info"][i][0].as<int>(),
                                     node["location_info"][i][1].as<float>(),
                                     node["location_info"][i][2].as<float>(), 0);
        }

        for (std::size_t loc = 0; loc < number_of_locations; loc++) {
            auto input_loc = node["age_distribution_by_location"].size() < number_of_locations ? 0 : loc;

            for (std::size_t i = 0; i < node["age_distribution_by_location"][input_loc].size(); i++) {
                location_db[loc].age_distribution.push_back(
                        node["age_distribution_by_location"][input_loc][i].as<double>());
            }
        }
        for (std::size_t loc = 0; loc < number_of_locations; loc++) {
            auto input_loc = node["p_treatment_for_less_than_5_by_location"].size() < number_of_locations ? 0 : loc;
            location_db[loc].p_treatment_less_than_5 = node["p_treatment_for_less_than_5_by_location"][input_loc].as<float>();
        }
        for (std::size_t loc = 0; loc < number_of_locations; loc++) {
            auto input_loc = node["p_treatment_for_more_than_5_by_location"].size() < number_of_locations ? 0 : loc;
            location_db[loc].p_treatment_more_than_5 = node["p_treatment_for_more_than_5_by_location"][input_loc].as<float>();
        }

        // If a raster was loaded for these items then use that instead
        if (SpatialData::get_instance().has_raster(SpatialData::SpatialFileType::Beta)) {
            printf("Beta raster and value supplied, ignoring beta_by_location setting\n");
        } else {
            for (std::size_t loc = 0; loc < number_of_locations; loc++) {
                auto input_loc = node["beta_by_location"].size() < number_of_locations ? 0 : loc;
                location_db[loc].beta = node["beta_by_location"][input_loc].as<float>();
            }
        }
        if (SpatialData::get_instance().has_raster(SpatialData::SpatialFileType::Population)) {
            printf("Population raster and value supplied, ignoring population_size_by_location setting\n");
        } else {
            for (std::size_t loc = 0; loc < number_of_locations; loc++) {
                auto input_loc = node["population_size_by_location"].size() < number_of_locations ? 0 : loc;
                location_db[loc].population_size = node["population_size_by_location"][input_loc].as<int>();
            }
        }

        return true;
    }
};

template <>
struct convert<ParasiteDensityLevel> {
  static Node encode(const ParasiteDensityLevel& rhs) {
    Node node;
    node.push_back("ParasiteDensityLevel");
    return node;
  }

  static bool decode(const Node& node, ParasiteDensityLevel& parasite_density_level) {
    //
    // if (!node.IsScalar()) {
    //   return false;
    parasite_density_level.log_parasite_density_cured = node["log_parasite_density_cured"].as<double>();
    parasite_density_level.log_parasite_density_from_liver = node["log_parasite_density_from_liver"].as<double>();
    parasite_density_level.log_parasite_density_asymptomatic = node["log_parasite_density_asymptomatic"].as<double>();
    parasite_density_level.log_parasite_density_clinical = node["log_parasite_density_clinical"].as<double>();
    parasite_density_level.log_parasite_density_clinical_from = node["log_parasite_density_clinical_from"].as<double>();
    parasite_density_level.log_parasite_density_clinical_to = node["log_parasite_density_clinical_to"].as<double>();
    parasite_density_level.log_parasite_density_detectable = node["log_parasite_density_detectable"].as<double>();

    parasite_density_level.log_parasite_density_detectable_pfpr =
        node["log_parasite_density_detectable_pfpr"] ? node["log_parasite_density_detectable_pfpr"].as<double>()
                                                     : node["log_parasite_density_detectable"].as<double>();
    parasite_density_level.log_parasite_density_pyrogenic = node["log_parasite_density_pyrogenic"].as<double>();

    return true;
  }
};

template <>
struct convert<RelativeInfectivity> {
  static Node encode(const RelativeInfectivity& rhs) {
    Node node;
    node.push_back("RelativeInfectivity");
    return node;
  }

  static bool decode(const Node& node, RelativeInfectivity& relative_infectivity) {
    relative_infectivity.sigma = node["sigma"].as<double>();
    const auto ro = node["ro"].as<double>();
    const auto blood_meal_volume = node["blood_meal_volume"].as<double>();

    const auto d_star = 1 / blood_meal_volume;

    relative_infectivity.ro_star = (log(ro) - log(d_star)) / relative_infectivity.sigma;

    relative_infectivity.sigma = log(10) / relative_infectivity.sigma;
    return true;
  }
};

template <>
struct convert<PfGeneInfo> {
  static Node encode(const PfGeneInfo& rhs) {
    Node node;
    node.push_back("GeneInfo");
    return node;
  }

  static bool decode(const Node& node, PfGeneInfo& gene_info) {
    for (const auto& chromosome_node : node) {
      auto chromosome = chromosome_node["chromosome"].as<int>();
      for (const auto& gene_node : chromosome_node["genes"]) {
        auto gene = gene_node.as<GeneInfo>();
        gene.chromosome = chromosome;
        gene_info.chromosome_infos[chromosome - 1].gene_infos.push_back(gene);
      }
    }
    return true;
  }
};

template <>
struct convert<GeneInfo> {
  static Node encode(const GeneInfo& rhs) {
    Node node;
    node.push_back("Gene");
    return node;
  }

  static bool decode(const Node& node, GeneInfo& gene) {
    gene.aa_position_infos.clear();
    gene.name = node["name"].as<std::string>();

    gene.max_copies = node["max_copies"] ? node["max_copies"].as<int>() : 1;
    gene.cnv_daily_crs =
        node["cnv_daily_crs"] ? node["cnv_daily_crs"].as<std::vector<double>>() : std::vector<double>();

    if (node["cnv_multiplicative_effect_on_EC50"]) {
      for (const auto& drug_node : node["cnv_multiplicative_effect_on_EC50"]) {
        gene.cnv_multiplicative_effect_on_EC50[drug_node["drug_id"].as<int>()] =
            drug_node["factors"].as<std::vector<double>>();
      }
    }

    if (node["multiplicative_effect_on_EC50_for_2_or_more_mutations"]) {
      for (const auto& drug_node : node["multiplicative_effect_on_EC50_for_2_or_more_mutations"]) {
        gene.multiplicative_effect_on_EC50_for_2_or_more_mutations[drug_node["drug_id"].as<int>()] =
            drug_node["factor"].as<double>();
      }
    }

    for (const auto& aa_node : node["aa_positions"]) {
      auto aa_position = aa_node.as<AaPositionInfo>();
      gene.aa_position_infos.push_back(aa_position);
    }
    return true;
  }
};

template <>
struct convert<AaPositionInfo> {
  static Node encode(const AaPositionInfo& rhs) {
    Node node;
    node.push_back("AaPosition");
    return node;
  }

  static bool decode(const Node& node, AaPositionInfo& aa_pos) {
    aa_pos.amino_acids.clear();
    aa_pos.daily_crs.clear();
    aa_pos.position = node["position"].as<int>();
    aa_pos.amino_acids = node["amino_acids"].as<std::vector<char>>();
    aa_pos.daily_crs = node["daily_crs"].as<std::vector<double>>();

    for (const auto& drug_node : node["multiplicative_effect_on_EC50"]) {
      aa_pos.multiplicative_effect_on_EC50[drug_node["drug_id"].as<int>()] =
          drug_node["factors"].as<std::vector<double>>();
    }

    return true;
  }
};

template <>
struct convert<OverrideEC50Patterns> {
  static Node encode(const OverrideEC50Patterns& rhs) {
    Node node;
    node.push_back("OverrideEC50Patterns");
    return node;
  }

  static bool decode(const Node& node, OverrideEC50Patterns& patterns) {
    for (const auto& pattern_node : node) {
      patterns.push_back(OverrideEC50Pattern { pattern_node["pattern"].as<std::string>(),
                                               pattern_node["drug_id"].as<int>(), pattern_node["ec50"].as<double>() });
    }
    return true;
  }
};

template <>
struct convert<MosquitoConfig> {
  static Node encode(const MosquitoConfig& mdb) {
    Node node;
    node.push_back("mosquito_config");
    return node;
  }
  static bool decode(const Node& node, MosquitoConfig& mcf) {
    mcf.prmc_size = node["prmc_size"].as<int>();
    if (node["interrupted_feeding_rate"]) {
      mcf.interrupted_feeding_rate = node["interrupted_feeding_rate"].as<std::vector<double>>();
    }
    if (node["interrupted_feeding_rate_raster"]) {
      mcf.interrupted_feeding_rate_raster = node["interrupted_feeding_rate"].as<std::string>();
    }
    if (!node["interrupted_feeding_rate"] && !node["interrupted_feeding_rate_raster"]) {
      LOG(FATAL) << "Either interrupted feeding rate or raster file needs to be supplied";
    }
    return true;
  }
};

}  // namespace YAML
#endif  // YAMLCONVERTER_H
