/* 
 * File:   Reporter.h
 * Author: Merlin
 *
 * Created on August 1, 2013, 12:05 PM
 */

#ifndef REPORTER_CUH
#define REPORTER_CUH

#include "Core/PropertyMacro.h"
#include <map>

class Model;

namespace GPU {
  class Reporter;
}

class GPU::Reporter {
 DISALLOW_COPY_AND_ASSIGN(Reporter)

 DISALLOW_MOVE(Reporter)

 POINTER_PROPERTY(Model, model)

  enum ReportType {
    CONSOLE,
    GUI,
    MONTHLY_REPORTER,
    MMC_REPORTER,
    TACT_REPORTER,
    NOVEL_DRUG_REPOTER,
    VALIDATION_REPORTER
  };
 public:
  static std::map<std::string, GPU::Reporter::ReportType> ReportTypeMap;

 public:
  Reporter();

  //    Reporter(const Reporter& orig);
  virtual ~Reporter();

  virtual void initialize() = 0;

  virtual void before_run() = 0;

  virtual void after_run() = 0;

  virtual void begin_time_step() = 0;

  virtual void monthly_report() = 0;

  static GPU::Reporter *MakeReport(GPU::Reporter::ReportType report_type);

 private:

};

#endif /* REPORTER_CUH */