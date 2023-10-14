/* 
 * File:   PkPdReporter.h
 * Author: Merlin
 *
 * Created on October 29, 2014, 12:56 PM
 */

#ifndef PKPDREPORTER_H
#define    PKPDREPORTER_H

#include "Core/PropertyMacro.h"
#include "Core/TypeDef.h"
#include "Reporters/Reporter.h"
#include <sstream>
#include <fstream>

class AppInput;

class PkPdReporter : public Reporter {
 DISALLOW_COPY_AND_ASSIGN(PkPdReporter)

 PROPERTY_REF(DoubleVector, yesterday_density)

 public:
    std::stringstream ss;
    const std::string group_sep = "-1111\t";
    const std::string sep = ",";

    PkPdReporter(AppInput* appInput=nullptr);

  virtual ~PkPdReporter();

  void initialize() override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  virtual void after_time_step();

  void monthly_report() override;

 private:
    AppInput* appInput{nullptr};
    std::ofstream outputFStream;

};

#endif    /* PKPDREPORTER_H */
