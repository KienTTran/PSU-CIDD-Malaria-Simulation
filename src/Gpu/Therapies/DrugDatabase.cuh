/* 
 * File:   DrugDatabase.h
 * Author: nguyentran
 *
 * Created on June 3, 2013, 11:05 AM
 */

#ifndef DRUGDATABASE_CUH
#define    DRUGDATABASE_CUH

#include "Core/PropertyMacro.h"
#include "DrugType.cuh"
#include <map>

namespace GPU{
  class DrugDatabase;
  typedef std::map<int, GPU::DrugType *> DrugTypePtrMap;
}



class GPU::DrugDatabase : public GPU::DrugTypePtrMap {
 DISALLOW_COPY_AND_ASSIGN(DrugDatabase)
  // VIRTUAL_PROPERTY_REF(DrugTypePtrMap, drug_db)

 public:
  DrugDatabase();

  //    DrugDatabase(const DrugDatabase& orig);
  virtual ~DrugDatabase();

  void add(GPU::DrugType *dt);


 private:

};

#endif    /* DRUGDATABASE_CUH */

