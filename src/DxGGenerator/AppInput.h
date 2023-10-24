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
};
