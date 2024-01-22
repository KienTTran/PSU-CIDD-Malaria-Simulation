/* 
 * File:   Reporter.cpp
 * Author: Merlin
 * 
 * Created on August 1, 2013, 12:05 PM
 */

#include "Reporter.cuh"
#include "ConsoleReporter.cuh"
#include "Model.h"
#include "MonthlyReporter.cuh"
#include "MMCReporter.cuh"
#include "TACTReporter.cuh"
#include "NovelDrugReporter.cuh"
#include "ValidationReporter.cuh"

std::map<std::string, GPU::Reporter::ReportType> GPU::Reporter::ReportTypeMap{
    {"Console",         CONSOLE},
    {"MonthlyReporter", MONTHLY_REPORTER},
    {"MMC",             MMC_REPORTER},
    {"TACT",            TACT_REPORTER},
    {"NovelDrug",       NOVEL_DRUG_REPOTER},
    {"ValidationReporter",       VALIDATION_REPORTER},
};

GPU::Reporter::Reporter() : model_(nullptr) {
}

GPU::Reporter::~Reporter() = default;

GPU::Reporter* GPU::Reporter::MakeReport(GPU::Reporter::ReportType report_type) {
  switch (report_type) {
    case CONSOLE:
      return new GPU::ConsoleReporter();
    case MONTHLY_REPORTER:
      return new GPU::MonthlyReporter();
    case MMC_REPORTER:
      return new GPU::MMCReporter();
    case TACT_REPORTER:
      return new GPU::TACTReporter();
    case NOVEL_DRUG_REPOTER:
      return new GPU::NovelDrugReporter();
    case VALIDATION_REPORTER:
      return new GPU::ValidationReporter();
    default:
      return new GPU::MonthlyReporter();
  }
}