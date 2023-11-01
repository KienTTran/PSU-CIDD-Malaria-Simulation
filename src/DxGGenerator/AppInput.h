#include <string>
#include <vector>

struct AppInput {
    std::string input_file { "input.yml" };
    std::string output_file { "" };
    std::vector<int> therapies{};
    std::vector<int> therapy_list{};
    std::vector<std::string> genotypes{};
    bool is_crt_calibration = false;
    double as_iov = -1.0;
    double as_iiv = -1.0;
    double as_ec50 = -1.0;
    bool is_ee_calibration = false;
    int number_of_drugs_in_combination { 1 };
    bool is_art { false };
    std::vector<int> dosing_days;
    std::vector<double> half_life;
    std::vector<double> k_max;
    std::vector<double> EC50;
    std::vector<double> slope;
    std::vector<double> mean_drug_absorption;
    int population_size { 10000 };
};
