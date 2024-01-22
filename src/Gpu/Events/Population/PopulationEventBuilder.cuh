#ifndef POPULATIONEVENTBUILDER_CUH
#define POPULATIONEVENTBUILDER_CUH

#include "Gpu/Events/Event.cuh"
#include <vector>

class Config;

namespace GPU{
    class PopulationEventBuilder;
}

namespace YAML {
class Node;
}

class GPU::PopulationEventBuilder {
public:
  static std::vector<GPU::Event*> build_introduce_parasite_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_parasites_periodically_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_parasites_periodically_events_v2(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_change_treatment_coverage_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_change_treatment_strategy_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_single_round_mda_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_modify_nested_mft_strategy_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_turn_on_mutation_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_turn_off_mutation_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_change_interrupted_feeding_rate_event(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_plas2_parasite_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_580Y_mutant_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_aq_mutant_parasite_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*>
  build_introduce_lumefantrine_mutant_parasite_events(const YAML::Node& node, Config* config);

  static std::vector<GPU::Event*> build_introduce_triple_mutant_to_dpm_parasite_events(
      const YAML::Node& node, Config* config
  );
  static std::vector<GPU::Event*> build_change_within_host_induced_free_recombination_events(const YAML::Node node, Config* pConfig);

  static std::vector<GPU::Event*> build_change_mutation_probability_by_locus_events(const YAML::Node node, Config* pConfig);
};

#endif // POPULATIONEVENTBUILDER_H
