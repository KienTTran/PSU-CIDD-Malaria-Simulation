/* 
 * File:   PersonIndexByLocationStateAgeClass.cu
 * Author: nguyentran
 * 
 * Created on May 2, 2013, 10:41 AM
 */

#include "PersonIndexByLocationStateAgeClass.cuh"
#include "PersonIndexByLocationStateAgeClassHandler.cuh"
#include "Core/Config/Config.h"
#include "Model.h"

#include <cassert>

GPU::PersonIndexByLocationStateAgeClass::PersonIndexByLocationStateAgeClass(const int &no_location, const int &no_host_state,
                                                                       const int &no_age_class) {
  Initialize(no_location, no_host_state, no_age_class);
}

GPU::PersonIndexByLocationStateAgeClass::~PersonIndexByLocationStateAgeClass() {

}

void GPU::PersonIndexByLocationStateAgeClass::Initialize(const int &no_location, const int &no_host_state,
                                                    const int &no_age_class) {
  vPerson_.clear();

  GPUPersonPtrVector ppv;

  GPUPersonPtrVector2 ppv2;
  ppv2.assign(no_age_class, ppv);

  GPUPersonPtrVector3 ppv3;
  ppv3.assign(no_host_state, ppv2);

  vPerson_.assign(no_location, ppv3);
}

void GPU::PersonIndexByLocationStateAgeClass::add(GPU::Person *p) {
  assert(vPerson_.size() > p->location() && p->location() >= 0);
  assert(vPerson_[p->location()].size() > p->host_state());
  assert(vPerson_[p->location()][p->host_state()].size() > p->age_class());
  assert(p->age_class() >= 0);

  add(p, p->location(), p->host_state(), p->age_class());

}

void GPU::PersonIndexByLocationStateAgeClass::add(GPU::Person *p, const int &location, const GPU::Person::HostStates &host_state,
                                             const int &age_class) {
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC add before loc {}->{} hs {}->{} ac {}->{} size {}->{}",
//                   p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size(),
//                   vPerson()[location][host_state][age_class].size());
  vPerson_[location][host_state][age_class].push_back(p);
//  printf("GPU::PersonIndexByLocationStateAgeClass::add loc %d hs %d ac %d person_index %d -> (size -1) %d\n",
//         location,host_state,age_class,
//         p->GPU::PersonIndexByLocationStateAgeClassHandler::index(),
//         vPerson_[location][host_state][age_class].size() - 1);
  p->GPU::PersonIndexByLocationStateAgeClassHandler::set_index(vPerson_[location][host_state][age_class].size() - 1);
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC add before loc {}->{} hs {}->{} ac {}->{} size {}->{}",
//                   p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size(),
//                   vPerson()[location][host_state][age_class].size());
}

void GPU::PersonIndexByLocationStateAgeClass::remove(GPU::Person *p) {
  remove_without_set_index(p);
  p->GPU::PersonIndexByLocationStateAgeClassHandler::set_index(-1);
}

void GPU::PersonIndexByLocationStateAgeClass::remove_without_set_index(GPU::Person *p) {
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC remove before loc {} hs {} ac {} size {}",
//                   p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size());
  vPerson_[p->location()][p->host_state()][p->age_class()].back()->GPU::PersonIndexByLocationStateAgeClassHandler::set_index(
      p->GPU::PersonIndexByLocationStateAgeClassHandler::index());
  vPerson_[p->location()][p->host_state()][p->age_class()][p->GPU::PersonIndexByLocationStateAgeClassHandler::index()] =
      vPerson_[p->location()][p->host_state()][p->age_class()].back();
  vPerson_[p->location()][p->host_state()][p->age_class()].pop_back();
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC remove after loc {} hs {} ac {} size {}",
//                   p->location(),static_cast<int>(p->host_state()),p->age_class(),
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size());
}

std::size_t GPU::PersonIndexByLocationStateAgeClass::size() const {
  return 0;
}

void
GPU::PersonIndexByLocationStateAgeClass::notify_change(GPU::Person *p, const GPU::Person::Property &property, const void *oldValue,
                                                  const void *newValue) {

  switch (property) {
    case GPU::Person::LOCATION: {
//      printf("GPU::PersonIndexByLocationStateAgeClass::notify_change LOCATION\n");
      change_property(p, *(int *) newValue, p->host_state(), p->age_class());
      break;
    }
    case GPU::Person::HOST_STATE: {
//      printf("GPU::PersonIndexByLocationStateAgeClass::notify_change HOST_STATE\n");
      change_property(p, p->location(), *(GPU::Person::HostStates *) newValue, p->age_class());
      break;
    }
    case GPU::Person::AGE_CLASS: {
//      printf("GPU::PersonIndexByLocationStateAgeClass::notify_change AGE_CLASS\n");
      change_property(p, p->location(), p->host_state(), *(int *) newValue);
      break;
    }
    default:break;
  }

}

void GPU::PersonIndexByLocationStateAgeClass::change_property(GPU::Person *p, const int &location,
                                                         const GPU::Person::HostStates &host_state, const int &age_class) {
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//  << "Changing " << p->location() << "-" << static_cast<int>(p->host_state()) << "-" << p->age_class() << " to "
//    << location << "-" << static_cast<int>(host_state) << "-" << age_class;
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC change_property before remove-add loc {}->{} hs {}->{} ac {}->{} size {}->{}",
//                   p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size(),
//                   vPerson()[location][host_state][age_class].size());
  //remove from old position
//  printf("GPU::PersonIndexByLocationStateAgeClass::change_property before remove p->loc %d -> %d, p->hs %d -> %d, p->ac %d -> %d\n",
//         p->location(),location,p->host_state(),host_state,p->age_class(),age_class);
//  printf("GPU::PersonIndexByLocationStateAgeClass::change_property before remove p->loc %d p->hs %d p->ac %d size %d person_index %d\n",
//         p->location(),p->host_state(),p->age_class(),
//         vPerson_[p->location()][p->host_state()][p->age_class()].size(),
//         p->GPU::PersonIndexByLocationStateAgeClassHandler::index());
  remove_without_set_index(p); //to save 1 set and improve performance since the index of p will changed when add
//  printf("GPU::PersonIndexByLocationStateAgeClass::change_property after remove p->loc %d p->hs %d p->ac %d size %d\n",
//         p->location(),p->host_state(),p->age_class(),
//         vPerson_[p->location()][p->host_state()][p->age_class()].size());
  //add to new position
  add(p, location, host_state, age_class);
//  LOG_IF(p->location() == 0 && p->host_state() == 2 && p->age_class() == 14,INFO)
//    << fmt::format("LocHostAC change_property after remove-add loc {}->{} hs {}->{} ac {}->{} size {}->{}",
//                   p->location(),location,static_cast<int>(p->host_state()),static_cast<int>(host_state),p->age_class(),age_class,
//                   vPerson()[p->location()][p->host_state()][p->age_class()].size(),
//                   vPerson()[location][host_state][age_class].size());
}

void GPU::PersonIndexByLocationStateAgeClass::update() {
  for (int location = 0; location < Model::CONFIG->number_of_locations(); location++) {
    for (int hs = 0; hs < GPU::Person::NUMBER_OF_STATE; hs++) {
      for (int ac = 0; ac < Model::CONFIG->number_of_age_classes(); ac++) {
        std::vector<GPU::Person *>(vPerson_[location][hs][ac]).swap(vPerson_[location][hs][ac]);
        vPerson_[location][hs][ac].shrink_to_fit();
      }
    }
  }

}
