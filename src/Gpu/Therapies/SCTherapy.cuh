/* 
 * File:   Therapy.h
 * Author: nguyentran
 *
 * Created on June 3, 2013, 7:50 PM
 */

#ifndef SCTHERAPY_CUH
#define    SCTHERAPY_CUH

#include "Core/PropertyMacro.h"
#include "Therapy.cuh"
#include <vector>

namespace GPU{
    class Therapy;
    class SCTherapy;
    class DrugType;
}

class GPU::SCTherapy : public GPU::Therapy {
DISALLOW_COPY_AND_ASSIGN(SCTherapy)

public:
  int artemisinin_id;
  std::vector<int> dosing_day;

public:
  SCTherapy();

  //    Therapy(const Therapy& orig);
  virtual ~SCTherapy();

  void add_drug(int drug_id);

  int get_arteminsinin_id() const;

  int get_max_dosing_day() const;
  //    int get_therapy_duration(int dosing_day);

private:

};

#endif    /* THERAPY_CUH */

