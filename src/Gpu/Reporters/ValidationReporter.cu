//
// Created by kient on 7/7/2022.
//

#include "ValidationReporter.cuh"
#include <date/date.h>
#include <iomanip>

#include "Constants.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Model.h"
#include "Gpu/Population/Population.cuh"
#include "Gpu/Population/Properties/PersonIndexByLocationStateAgeClass.cuh"
#include "Gpu/Strategies/IStrategy.cuh"
#include "ReporterUtils.cuh"
#include "easylogging++.h"

GPU::ValidationReporter::ValidationReporter() = default;

GPU::ValidationReporter::~ValidationReporter() = default;

void GPU::ValidationReporter::initialize() {
    monthly_data_file.open(fmt::format("{}/validation_monthly_data_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
    summary_data_file.open(fmt::format("{}/validation_summary_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
    gene_freq_file.open(fmt::format("{}/validation_gene_freq_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
    gene_db_file.open(fmt::format("{}/validation_gene_db_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
    prmc_freq_file.open(fmt::format("{}/validation_prmc_freq_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
    prmc_db_file.open(fmt::format("{}/validation_prmc_db_{}.txt", Model::MODEL->output_path(), Model::MODEL->cluster_job_number()));
}

void GPU::ValidationReporter::before_run() {}

void GPU::ValidationReporter::begin_time_step() {}

void GPU::ValidationReporter::monthly_report() {
    std::stringstream ss;

    ss << Model::GPU_SCHEDULER->current_time() << sep;
    ss << std::chrono::system_clock::to_time_t(Model::GPU_SCHEDULER->calendar_date) << sep;
    ss << date::format("%Y\t%m\t%d", Model::GPU_SCHEDULER->calendar_date) << sep;
    ss << Model::CONFIG->seasonal_info()->get_seasonal_factor(Model::GPU_SCHEDULER->calendar_date, 0) << sep;
    ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 1) << sep;
    ss << Model::TREATMENT_COVERAGE->get_probability_to_be_treated(0, 10) << sep;
    ss << Model::GPU_POPULATION->size() << sep;
    ss << group_sep;//10
    print_EIR_PfPR_by_location(ss);//12
    ss << group_sep;//16
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_new_infections_by_location()[loc] << sep;
    }
    ss << group_sep;//18
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_treatment_by_location()[loc] << sep; //Incidence
    }
    ss << group_sep;//20
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] << sep;
    }
    ss << group_sep;///22
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for(auto moi : Model::GPU_DATA_COLLECTOR->multiple_of_infection_by_location()[loc]){
            ss << moi << sep;
        }
    }
    ss << group_sep;///33
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (auto ac = 0ul; ac < Model::CONFIG->number_of_age_classes(); ac++){
            ss << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location_age_group()[loc][ac] << sep;
        }
    }
    ss << group_sep;//49
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (int age = 0; age < 80; age++){
            ss << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location_age()[loc][age] << sep;
        }
    }
    ss << group_sep;///130
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->current_TF_by_location()[loc] << sep;
    }
    ss << group_sep;//132
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_mutation_events_by_location()[loc] << sep;
    }
    ss << group_sep;//134
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (auto ac = 0ul; ac < Model::CONFIG->number_of_age_classes(); ac++){
            ss << Model::GPU_DATA_COLLECTOR->number_of_treatments_by_location_age_year()[loc][ac] << sep;
        }
    }
    ss << group_sep;//150
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->total_number_of_bites_by_location()[loc] << sep;
    }
    ss << group_sep;//152
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->total_number_of_bites_by_location_year()[loc] << sep;
    }
    ss << group_sep;//154
    ss << Model::GPU_DATA_COLLECTOR->number_of_treatments_with_therapy_ID()[0] << sep;
    ss << group_sep;//156
    ss << Model::GPU_DATA_COLLECTOR->number_of_treatments_fail_with_therapy_ID()[0] << sep;
    ss << group_sep;//158
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        if ((Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::ASYMPTOMATIC] + Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::CLINICAL]) == 0){
            ss << 0 << sep;
        }
        else{
            ss << Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::ASYMPTOMATIC] /
                  (Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::ASYMPTOMATIC] + Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::CLINICAL])  << sep;
        }
    }
    ss << group_sep;//160
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (auto age = 0; age < 80; age++){
            ss << Model::GPU_DATA_COLLECTOR->popsize_by_location_age()[loc][age] << sep;
        }
    }
    ss << group_sep;//241
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (auto age = 0; age < 80; age++){
            ss << Model::GPU_DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location_age()[loc][age] << sep;
        }
    }
    ss << group_sep;//322
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        for (auto ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++){
            int all_infected_pop = Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate_age_class()[loc][GPU::Person::ASYMPTOMATIC][ac]
                                 + Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate_age_class()[loc][GPU::Person::CLINICAL][ac];
            if (all_infected_pop == 0){
                ss << 0 << sep;
            }
            else{
                ss << std::setprecision(6) << Model::GPU_DATA_COLLECTOR->number_of_clinical_by_location_age_group()[loc][ac] / (double) all_infected_pop<< sep;
            }
        }
    }
    ss << group_sep;//338
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->total_immune_by_location()[loc] / Model::GPU_POPULATION->size(loc) << sep;
    }
    ss << group_sep;//340
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        int all_infected_pop = Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::ASYMPTOMATIC]
                               + Model::GPU_DATA_COLLECTOR->popsize_by_location_hoststate()[loc][GPU::Person::CLINICAL];
        if (all_infected_pop == 0){
            ss << 0 << sep;
        }
        else{
            ss << std::setprecision(6) << Model::GPU_DATA_COLLECTOR->monthly_number_of_clinical_episode_by_location()[loc] / (double) all_infected_pop<< sep;
        }
    }
    ss << group_sep;//342
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        ss << Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[loc] << sep;
    }
    ss << group_sep;

    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); loc++) {
        LOG_IF(Model::CONFIG->debug_config().enable_debug_text, INFO) << fmt::format("ValidationReporter::monthly_report: {} PfPR {}",
                                                                                     loc,
                                                                                     Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[loc] * 100.0f);
    }
    LOG_IF(Model::CONFIG->debug_config().enable_debug_text, INFO) << fmt::format("ValidationReporter::monthly_report: TF {} {} {}",
                                                                                 Model::GPU_DATA_COLLECTOR->current_TF_by_location()[6],
                                                                                 Model::GPU_DATA_COLLECTOR->current_TF_by_location()[7],
                                                                                 Model::GPU_DATA_COLLECTOR->current_TF_by_location()[8]);

    monthly_data_file << ss.str() << std::endl;

    std::stringstream gene_freq_ss;
//    GPU::ReporterUtils::output_genotype_frequency3(gene_freq_ss, Model::CONFIG->gpu_genotype_db.size(),
//                                              Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>());

    std::stringstream prmc_freq_ss;
    GPU::ReporterUtils::output_genotype_frequency4(gene_freq_ss, prmc_freq_ss, Model::CONFIG->gpu_genotype_db.size(),
                                              Model::GPU_POPULATION->get_person_index<GPU::PersonIndexByLocationStateAgeClass>());

    gene_freq_file << gene_freq_ss.str() << std::endl;
    prmc_freq_file << prmc_freq_ss.str() << std::endl;
}

void GPU::ValidationReporter::after_run() {
    std::stringstream ss;

    ss.str("");
    ss << Model::RANDOM->seed() << sep << Model::CONFIG->number_of_locations() << sep;
    ss << Model::CONFIG->location_db()[0].beta << sep;
    ss << Model::CONFIG->location_db()[0].population_size << sep;
    print_EIR_PfPR_by_location(ss);
    ss << group_sep;//10
    // output last strategy information
    ss << Model::GPU_TREATMENT_STRATEGY->id << sep;//11
    // output NTF
    const auto total_time_in_years = (Model::GPU_SCHEDULER->current_time() - Model::CONFIG->start_of_comparison_period())
                                     / static_cast<double>(Constants::DAYS_IN_YEAR());
    auto sum_ntf = 0.0;
    ul pop_size = 0;
    for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++) {
        sum_ntf += Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[location];
        pop_size += Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];
    }
    ss << (sum_ntf * 100 / pop_size) / total_time_in_years << sep;
    ss << group_sep;//13
    for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
        ss << Model::GPU_DATA_COLLECTOR->cumulative_clinical_episodes_by_location()[location] << sep;
//        std::cout << Model::GPU_DATA_COLLECTOR->cumulative_clinical_episodes_by_location()[location] << "\t";
    }
    ss << group_sep;//15
    for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
        ss << Model::GPU_DATA_COLLECTOR->percentage_bites_on_top_20_by_location()[location] * 100 << "%"  << sep;
//        std::cout << Model::GPU_DATA_COLLECTOR->percentage_bites_on_top_20_by_location()[location] * 100 << "%" << "\t";
    }
    ss << group_sep;//17
    for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
        double location_NTF = Model::GPU_DATA_COLLECTOR->cumulative_NTF_by_location()[location] * 100 /
                              (double) Model::GPU_DATA_COLLECTOR->popsize_by_location()[location];
        location_NTF /= total_time_in_years;
        ss << location_NTF << sep;
//        std::cout << location_NTF << "\t";
    }
    ss << group_sep;//19
    for (auto location = 0; location < Model::CONFIG->number_of_locations(); location++)
    {
        for (auto age = 0; age < 80; age++){
            ss << Model::GPU_DATA_COLLECTOR->cumulative_clinical_episodes_by_location_age()[location][age]/total_time_in_years/Model::GPU_DATA_COLLECTOR->popsize_by_location_age()[location][age] << sep;
        }
    }
    ss << group_sep;//100

    summary_data_file << ss.str() << std::endl;

    for (auto [g_id, genotype] : Model::CONFIG->gpu_genotype_db) {
        gene_db_file << g_id << sep << genotype->aa_sequence << std::endl;
        prmc_db_file << g_id << sep << genotype->aa_sequence << std::endl;
    }

    gene_db_file.close();
    gene_freq_file.close();
    prmc_db_file.close();
    prmc_freq_file.close();
    monthly_data_file.close();
    summary_data_file.close();
}

void GPU::ValidationReporter::print_EIR_PfPR_by_location(std::stringstream& ss) {
    for (auto loc = 0; loc < Model::CONFIG->number_of_locations(); ++loc) {
        //
        // EIR
        if (Model::GPU_DATA_COLLECTOR->EIR_by_location_year()[loc].empty()) {
            ss << 0 << sep;
        } else {
            ss << Model::GPU_DATA_COLLECTOR->EIR_by_location_year()[loc].back() << sep;
        }
        ss << group_sep;//11
        // pfpr <5 , 2-10 and all
        ss << Model::GPU_DATA_COLLECTOR->get_blood_slide_prevalence(loc, 2, 10) * 100 << sep;
        ss << Model::GPU_DATA_COLLECTOR->get_blood_slide_prevalence(loc, 0, 5) * 100 << sep;
        ss << Model::GPU_DATA_COLLECTOR->blood_slide_prevalence_by_location()[loc] * 100 << sep;
    }
}