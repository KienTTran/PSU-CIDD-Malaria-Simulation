/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.cpp
 * Author: Merlin
 *
 * Created on January 12, 2017, 4:31 PM
 */
// TODO: make it works, input for genotype will be string or a file with a list of string


#include <math.h>  // log10

#include <CLI/CLI.hpp>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "AppInput.h"
#include "Core/Config/Config.h"
#include "Core/Random.h"
#include "Events/ProgressToClinicalEvent.h"
#include "IndividualsFileReporter.h"
#include "MDC/ModelDataCollector.h"
#include "Model.h"
#include "PkPdReporter.h"
#include "Population/ImmuneSystem.h"
#include "Population/Person.h"
#include "Population/Population.h"
#include "Population/Properties/PersonIndexAll.h"
#include "Strategies/IStrategy.h"
#include "Strategies/SFTStrategy.h"
#include "Therapies//SCTherapy.h"
#include "Therapies/Therapy.h"
#include "Parasites/Genotype.h"
#include "Mosquito/Mosquito.h"
#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

using namespace std;

bool validate_config_for_ee(AppInput& input);
double getEfficacyForTherapy(std::string g_str, Model* p_model,AppInput& input, int therapy_id);
double getEfficacyForTherapyCRT(Model* p_model,AppInput& input, int therapy_id);

void create_cli_option(CLI::App& app,AppInput& input);

// efficacy_map efficacies;


inline double round(double val) {
    if (val < 0) return ceil(val - 0.5);
    return floor(val + 0.5);
}

/*
 *
 */
int main(int argc, char** argv) {
    CLI::App app { "PKPD model" };
    AppInput input;
    create_cli_option(app,input);
    CLI11_PARSE(app, argc, argv);

    // Turn off logger
    el::Configurations default_conf;
    default_conf.setToDefault();
    default_conf.setGlobally(el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", default_conf);
    START_EASYLOGGINGPP(argc, argv);

    auto p_model = new Model();
    p_model->set_config_filename(input.input_file);
    p_model->initialize();

    if (input.as_iov != -1) {
        p_model->CONFIG->as_iov() = input.as_iov;
    }

    if (input.as_iiv != -1) {
        for (auto& sd : p_model->CONFIG->drug_db()->at(0)->age_group_specific_drug_concentration_sd()) {
            sd = input.as_iiv;
        }
    }

    if (input.as_ec50 != -1) {
        // TODO: fix it
        //    p_model->CONFIG->EC50_power_n_table()[0][0] = pow(as_ec50, p_model->CONFIG->drug_db()->at(0)->n());
    }
    std::cout << std::setprecision(5);
    int max_therapy_id { 0 }, min_therapy_id { 0 };

    if(input.therapy_list.empty()){
        if (input.therapies.empty()) {
            min_therapy_id = 0;
            max_therapy_id = 0;
        } else if (input.therapies.size() == 1) {
            min_therapy_id = input.therapies[0];
            max_therapy_id = input.therapies[0];
        } else if (input.therapies.size() == 2) {
            min_therapy_id = input.therapies[0];
            max_therapy_id = input.therapies[1];
        }
    }

    // TODO: Genotype should be imported  from input files

    if(input.is_crt_calibration){
        if(input.genotypes.empty()){
            std::cout << "List of population genotypes is empty" << std::endl;
            exit(0);
        }
        std::stringstream ss;
        if(input.therapy_list.empty()){
            for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
              std::cout << *Model::CONFIG->therapy_db()[therapy_id] << "\t" ;
            }
        }
        else{
            for (auto therapy_id : input.therapy_list) {
              std::cout << *Model::CONFIG->therapy_db()[therapy_id] << "\t";
            }
        }
        std::cout << std::endl;
        if(input.therapy_list.empty()){
            for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
              double efficacy = getEfficacyForTherapyCRT(p_model, input, therapy_id);
              ss << efficacy << (therapy_id == max_therapy_id ? "" : "\t");
            }
        }
        else{
            for (int t_index = 0; t_index < input.therapy_list.size(); t_index++) {
              double efficacy = getEfficacyForTherapyCRT(p_model, input, input.therapy_list[t_index]);
              ss << efficacy << (input.therapy_list[t_index] == input.therapy_list.size() - 1 ? "" : "\t");
            }
        }
        std::cout << ss.str() << std::endl;
    }
    else if(input.is_ee_calibration){
        if (!validate_config_for_ee(input)){
            std::cout << "Parameters for Efficacy Estimator are not correct" << std::endl;
            exit(0);
        }
        else if(input.genotypes.size() > 1){
              std::cout << "Only 1 genotype is accepted using Efficacy Estimator" << std::endl;
              exit(0);
        }
        else{
            // ==== override population size ======
            if (input.population_size != Model::POPULATION->size()) {
                Model::CONFIG->location_db()[0].population_size = input.population_size;
                delete Model::POPULATION;
                p_model->set_population(new Population(p_model));
                Model::POPULATION = p_model->population();
                p_model->population()->initialize();
            }
            // ==== override drug type info ========
            auto start_drug_id = input.is_art ? 0 : 1;
            for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
              auto* dt = Model::CONFIG->drug_db()->at(i + start_drug_id);
              dt->set_name(fmt::format("D{}", i));
              dt->set_drug_half_life(input.half_life[i]);
              dt->set_maximum_parasite_killing_rate(input.k_max[i]);
              dt->set_n(input.slope[i]);
              // TODO: add app arguments later
              //    dt->set_p_mutation(0.0);
              dt->set_k(4);
              for (double& mda:dt->age_specific_drug_absorption()) {
                mda = input.mean_drug_absorption[i];
              }
              //    Model::CONFIG->EC50_power_n_table()[0][i + start_drug_id] = pow(input.EC50[i], dt->n());
            }

            // ======= override therapy 0 ==========
            auto scTherapy = dynamic_cast<SCTherapy*>(Model::CONFIG->therapy_db()[0]);
            scTherapy->drug_ids.clear();
            scTherapy->dosing_day.clear();

            for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
              scTherapy->drug_ids.push_back(i + start_drug_id);
              scTherapy->dosing_day.push_back(input.dosing_days[i]);
            }

            // ==========reset and override reporters ==================
            for (auto reporter : p_model->reporters()) {
              delete reporter;
            }
            p_model->reporters().clear();
            p_model->add_reporter(new PkPdReporter(&input));

            Model::CONFIG->genotype_db.clear();
            std::vector<Genotype*> genotype_inputs;
            for(auto genotype_str : input.genotypes){
              genotype_inputs.push_back(Model::CONFIG->genotype_db.get_genotype(genotype_str,p_model->CONFIG));
            }

            // =========infect population with genotype 0================
            auto* genotype = Model::CONFIG->genotype_db.get_genotype(genotype_inputs.front()->aa_sequence,p_model->CONFIG);

            for (auto person : Model::POPULATION->all_persons()->vPerson()) {
              auto density = Model::CONFIG->parasite_density_level().log_parasite_density_from_liver;
              auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

              person->immune_system()->set_increase(true);
              person->set_host_state(Person::EXPOSED);

              blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
              blood_parasite->set_last_update_log10_parasite_density(density);

              ProgressToClinicalEvent::schedule_event(Model::SCHEDULER, person, blood_parasite, 0);
            }

            // run model
            p_model->run();

            const auto result = 1 - Model::DATA_COLLECTOR->blood_slide_prevalence_by_location()[0];
            fmt::print(
                "pop\tdose\thalflife\tkmax\tec50\tslope\tis_art\tefficacy\n"
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:f}\n", p_model->population()->all_persons()->size(), fmt::join(input.dosing_days, "\t"),
                fmt::join(input.half_life, "\t"), fmt::join(input.k_max, "\t"), fmt::join(input.EC50, "\t"),
                fmt::join(input.slope, "\t"), input.is_art ? 1 : 0, result
            );
        }
    }
    else{
        std::cout << "ID\tGenotype\t";
        if(input.therapy_list.empty()){
            for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
              std::cout << *Model::CONFIG->therapy_db()[therapy_id] << "\t" ;
            }
        }
        else{
            for (auto therapy_id : input.therapy_list) {
              std::cout << *Model::CONFIG->therapy_db()[therapy_id] << "\t";
            }
        }
        std::cout << std::endl;
        for(int g_index = 0; g_index < input.genotypes.size(); g_index++){
            std::stringstream ss;
            ss << g_index << "\t" << Model::MOSQUITO->get_old_genotype_string2(input.genotypes[g_index]) << "\t";
            if(input.therapy_list.empty()){
              for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
                double efficacy = getEfficacyForTherapy(input.genotypes[g_index], p_model, input, therapy_id);
//                ss << efficacy << (therapy_id == max_therapy_id ? "" : "\t");
                  ss << efficacy << "\t";
              }
            }
            else{
              for (int t_index = 0; t_index < input.therapy_list.size(); t_index++) {
                double efficacy = getEfficacyForTherapy(input.genotypes[g_index], p_model, input, input.therapy_list[t_index]);
//                ss << efficacy << (input.therapy_list[t_index] == input.therapy_list.size() - 1 ? "" : "\t");
                ss << efficacy << "\t";
              }
            }
            std::cout << ss.str() << std::endl;
        }
    }
    delete p_model;
    return 0;
}

void create_cli_option(CLI::App& app, AppInput& input) {
    app.add_option("-i", input.input_file, "Input filename for DxG");
    app.add_option("-g", input.genotypes, "Genotype patterns for population (3 only) ex: [WT KEL1 KEL1/PLA1]\n"
                                          "Or genotype list without -cc flags");
    app.add_option("-t", input.therapies, "Get efficacy for range therapies [from to]");
    app.add_option("--of", input.output_file, "Output density to file");
    app.add_option("--iov", input.as_iov, "AS inter-occasion-variability");
    app.add_option("--iiv", input.as_iiv, "AS inter-individual-variability");
    app.add_option("--ec50", input.as_ec50, "EC50 for AS on C580 only");
    app.add_option("--pil", input.is_print_immunity_level, "print immunity level");
    //DxG to calibrate PfCRT factor. --cc to enable and use with -g to distribute 3 genotypes to population, --tl to get efficacy of list of therapy
    app.add_option("--cc", input.is_crt_calibration, "Enable PfCRT ec50 calibration");
    app.add_option("--tl", input.therapy_list, "Get efficacy for list of therapies [0 1 2 ...]");
    //EfficacyEstimator. --ee to enable and use below params
    app.add_option("--popsize", input.population_size, "Population size for EfficacyEstimator");
    app.add_option("--ee", input.is_ee_calibration, "Enable EfficacyEstimator");
    app.add_option("--EC50", input.EC50, "ee ec50");
    app.add_option("--art", input.is_art, "ee is art");
    app.add_option("--nd", input.number_of_drugs_in_combination, "ee number of drug");
    app.add_option("--dose", input.dosing_days, "ee dose");
    app.add_option("--kmax", input.k_max, "ee kmax");
    app.add_option("--halflife", input.half_life, "ee half life");
    app.add_option("--slope", input.slope, "ee slope");
    app.add_option("--mda", input.mean_drug_absorption, "ee mean drug absorption");
}

double getEfficacyForTherapy(std::string g_str, Model* p_model, AppInput& input, int therapy_id) {
    auto* mainTherapy = p_model->CONFIG->therapy_db()[therapy_id];
    dynamic_cast<SFTStrategy*>(Model::TREATMENT_STRATEGY)->get_therapy_list().clear();
    dynamic_cast<SFTStrategy*>(Model::TREATMENT_STRATEGY)->add_therapy(mainTherapy);

    // reset reporter
    for (auto reporter : p_model->reporters()) {
        delete reporter;
    }

    p_model->reporters().clear();

    if(!input.output_file.empty()){
        p_model->add_reporter(new PkPdReporter(&input));
    }
    else{
        p_model->add_reporter(new PkPdReporter());
    }

    for (auto person : Model::POPULATION->all_persons()->vPerson()) {
        auto density = Model::CONFIG->parasite_density_level().log_parasite_density_from_liver;
        auto* genotype = Model::CONFIG->genotype_db.get_genotype(g_str,p_model->CONFIG);
        auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

        person->immune_system()->set_increase(true);
        person->set_host_state(Person::EXPOSED);

        blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
        blood_parasite->set_last_update_log10_parasite_density(density);

        ProgressToClinicalEvent::schedule_event(Model::SCHEDULER, person, blood_parasite, 0);
    }

    p_model->run();
    const auto result = 1 - Model::DATA_COLLECTOR->blood_slide_prevalence_by_location()[0];

    delete Model::POPULATION;
    delete Model::SCHEDULER;
    p_model->set_population(new Population(p_model));
    Model::POPULATION = p_model->population();
    p_model->set_scheduler(new Scheduler(p_model));
    Model::SCHEDULER = p_model->scheduler();

    p_model->scheduler()->initialize(Model::CONFIG->starting_date(), Model::CONFIG->total_time());
    p_model->population()->initialize();

    return result;
}

double getEfficacyForTherapyCRT(Model* p_model, AppInput& input, int therapy_id) {
    auto* mainTherapy = Model::CONFIG->therapy_db()[therapy_id];
    dynamic_cast<SFTStrategy*>(Model::TREATMENT_STRATEGY)->get_therapy_list().clear();
    dynamic_cast<SFTStrategy*>(Model::TREATMENT_STRATEGY)->add_therapy(mainTherapy);

    // reset reporter
    for (auto reporter : p_model->reporters()) {
        delete reporter;
    }

    p_model->reporters().clear();

    if(!input.output_file.empty()){
        p_model->add_reporter(new PkPdReporter(&input));
    }
    else{
        p_model->add_reporter(new PkPdReporter());
    }

    for (auto person : Model::POPULATION->all_persons()->vPerson()) {
        //The genotype distribution is from table 2 in http://dx.doi.org/10.1016/S1473-3099(19)30391-3
        //Run these 3 genotypes in population with and without F145I,T93S,H97Y and I218F
        //to get ec50 of WT and mutant genotypes
        std::string g_str = "";
        int infect_prob = Model::RANDOM->random_uniform_int(1, 104);
        if(infect_prob < 74) {//KEL1/PLA1
            g_str = input.genotypes[2];
        }
        else if(infect_prob < 91){//KEL1
            g_str = input.genotypes[1];
        }
        else{//WT
            g_str = input.genotypes[0];
        }
        auto* genotype = Model::CONFIG->genotype_db.get_genotype(g_str,p_model->CONFIG);
        auto* blood_parasite = person->add_new_parasite_to_blood(genotype);
        auto density = Model::CONFIG->parasite_density_level().log_parasite_density_from_liver;

        person->immune_system()->set_increase(true);
        person->set_host_state(Person::EXPOSED);

        blood_parasite->set_gametocyte_level(Model::CONFIG->gametocyte_level_full());
        blood_parasite->set_last_update_log10_parasite_density(density);

        ProgressToClinicalEvent::schedule_event(Model::SCHEDULER, person, blood_parasite, 0);
    }

    p_model->run();
    const auto result = 1 - Model::DATA_COLLECTOR->blood_slide_prevalence_by_location()[0];

    delete Model::POPULATION;
    delete Model::SCHEDULER;
    p_model->set_population(new Population(p_model));
    Model::POPULATION = p_model->population();
    p_model->set_scheduler(new Scheduler(p_model));
    Model::SCHEDULER = p_model->scheduler();

    p_model->scheduler()->initialize(Model::CONFIG->starting_date(), Model::CONFIG->total_time());
    p_model->population()->initialize();

    return result;
}

bool validate_config_for_ee(AppInput& input) {
    input.number_of_drugs_in_combination = input.half_life.size();

    if (input.number_of_drugs_in_combination > 5) {
        std::cerr << "Error: Number of drugs in combination should not greater than 5" << std::endl;
        return false;
    }

    if (input.k_max.size() != input.number_of_drugs_in_combination
        || input.EC50.size() != input.number_of_drugs_in_combination
        || input.slope.size() != input.number_of_drugs_in_combination
        || input.dosing_days.size() != input.number_of_drugs_in_combination) {
        std::cerr << "Error: Wrong number of drugs in combination" << std::endl;
        return false;
    }

    for (auto k : input.k_max) {
        if (k >= 1 || k < 0) {
            std::cerr << "Error: k_max should be in range of (0,1]" << std::endl;
            return false;
        }
    }

    for (auto ec50 : input.EC50) {
        if (ec50 < 0) {
            std::cerr << "Error: EC50 should be greater than 0." << std::endl;
            return false;
        }
    }

    for (auto n : input.slope) {
        if (n < 0) {
            std::cerr << "Error: n should greater than 0." << std::endl;
            return false;
        }
    }

    for (auto dosing : input.dosing_days) {
        if (dosing < 0) {
            std::cerr << "Error: dosing should greater than 0." << std::endl;
            return false;
        }
    }

    if (input.mean_drug_absorption.empty()) {
        for (int i = 0; i < input.number_of_drugs_in_combination; ++i) {
            input.mean_drug_absorption.push_back(1.0);
        }
    } else if (input.mean_drug_absorption.size() != input.number_of_drugs_in_combination) {
        std::cerr << "Error: Wrong number of drugs in combination" << std::endl;
        return false;
    }

    return true;
}