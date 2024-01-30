/* 
 * File:   Therapy.h
 * Author: nguyentran
 *
 * Created on June 3, 2013, 7:50 PM
 */

#ifndef THERAPY_CUH
#define    THERAPY_CUH

#include "Core/PropertyMacro.h"
#include <vector>
#include <Model.h>
#include "DrugDatabase.cuh"
#include "Core/Config/Config.h"
#include "DrugType.cuh"

namespace GPU{
    class Therapy;
    class DrugType;
}

class GPU::Therapy {
    
public:
    DISALLOW_COPY_AND_ASSIGN(Therapy)

    VIRTUAL_PROPERTY_REF(int, id)

    VIRTUAL_PROPERTY_REF(int, testing_day)

    VIRTUAL_PROPERTY_REF(std::string, name)

public:
  std::vector<int> drug_ids;

public:
  Therapy();

  //    Therapy(const Therapy& orig);
  virtual ~Therapy();

  virtual void add_drug(int drug_id);

  friend std::ostream &operator<<(std::ostream &os, const Therapy &therapy){
      if ( !therapy.name_.empty() ){
        os << therapy.name_;
      } else {
        os << Model::CONFIG->gpu_drug_db()->at(therapy.drug_ids[0])->name();
        for (int i = 1; i < therapy.drug_ids.size(); ++i) {
          os << "+" << Model::CONFIG->gpu_drug_db()->at(therapy.drug_ids[i])->name();
        }
      }

      return os;
    }

private:

};

#endif    /* THERAPY_CUH */

