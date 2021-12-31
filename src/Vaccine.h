/*
  This file is part of the FRED system.

  Copyright (c) 2010-2015, University of Pittsburgh, John Grefenstette,
  Shawn Brown, Roni Rosenfield, Alona Fyshe, David Galloway, Nathan
  Stone, Jay DePasse, Anuroop Sriram, and Donald Burke.

  Licensed under the BSD 3-Clause license.  See the file "LICENSE" for
  more information.
*/

//
//
// File: VaccineDose.h
//

#ifndef _FRED_VACCINE_H
#define _FRED_VACCINE_H

#include "Global.h"
#include "Disease_List.h"
#include "Disease.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

class Vaccine_Dose;

class Vaccine{
public:
  // Creation
  Vaccine(string _name, int _id, int _disease, 
          int _total_avail, int _additional_per_day, 
          int _start_day, int num_strains, int *_strains);
  ~Vaccine();
  
  void add_dose(Vaccine_Dose* dose);
  void add_booster(Vaccine_Dose* dose);
  
  int get_disease()                 const { return disease; }
  int get_ID()                     const { return id; }
  int get_number_doses()           const { return doses.size(); }
  Vaccine_Dose* get_dose(int i)    const { return doses[i]; }
  int get_number_boosters()           const { return boosters.size(); }
  Vaccine_Dose* get_booster(int i)    const { return boosters[i]; }
  
  // Logistics Functions
  int get_initial_stock()          const { return initial_stock; }
  int get_total_avail()            const { return total_avail; }
  int get_current_reserve()        const { return reserve; }
  int get_current_stock()          const { return stock; }
  int get_additional_per_day()     const { return additional_per_day; }
  void add_stock( int add ){ 
    if(add <= reserve){
      stock   += add;
      reserve -= add;
    }
    else{
      stock   += reserve;
      reserve  = 0;
    }
  }
  
  void remove_stock( int remove ) {
    stock-=remove;
    if(stock < 0) stock = 0;
  }

  int get_num_strains() {
    return num_strains;
  }
  
  int get_strain(int i);
  
  void set_disease_specific_efficacy(int i, double eff_){
    if(i < disease_specific_efficacy_modifier.size() && eff_ <= 1.0 && eff_ >= 0){
      this->disease_specific_efficacy_modifier[i] = eff_;
    }
  }

  double get_disease_specific_efficacy(int i){
    if(i < disease_specific_efficacy_modifier.size()){
       return disease_specific_efficacy_modifier[i];
    }else{
      return -1;
    }
  }

  void set_disease_specific_efficacy_symp(int i, double eff_){
    if(i < disease_specific_efficacy_symp_modifier.size() && eff_ >= 0){
      this->disease_specific_efficacy_symp_modifier[i] = eff_;
    }
  }

  double get_disease_specific_efficacy_symp(int i){
    if(i < disease_specific_efficacy_symp_modifier.size()){
       return disease_specific_efficacy_symp_modifier[i];
    }else{
      return -1;
    }
  }

  void set_disease_specific_efficacy_hosp(int i, double eff_){
    if(i < disease_specific_efficacy_hosp_modifier.size() && eff_ >= 0){
      this->disease_specific_efficacy_hosp_modifier[i] = eff_;
    }
  }

  double get_disease_specific_efficacy_hosp(int i){
    if(i < disease_specific_efficacy_hosp_modifier.size()){
       return disease_specific_efficacy_hosp_modifier[i];
    }else{
      return -1;
    }
  }
  
  
  void set_disease_specific_efficacy(){
    //this->disease_specific_efficacy_modifier.clear();
    for(int dis_id = 0; dis_id < Global::Diseases.get_number_of_diseases(); ++dis_id){
      this->disease_specific_efficacy_modifier.push_back(1.0);
      this->disease_specific_efficacy_symp_modifier.push_back(1.0);
      this->disease_specific_efficacy_hosp_modifier.push_back(1.0);
    }
  }
  
  //Utility Functions
  void print() const;
  void update(int day);
  void update(int day, int add_vac);
  void reset();
  
private:
  string name;
  int id;                              // Which in the number of vaccines is it
  int disease;                          // Which Disease is this vaccine for
  int number_doses;                    // How many doses does the vaccine need.  
  vector < Vaccine_Dose* > doses;       // Data structure to hold the efficacy of each dose.
  int number_boosters;                    // How many boosters does the vaccine need.  
  vector < Vaccine_Dose* > boosters;       // Data structure to hold the efficacy of each booster.
  vector < double > disease_specific_efficacy_modifier;
  vector < double > disease_specific_efficacy_symp_modifier;
  vector < double > disease_specific_efficacy_hosp_modifier;
  
  
  int initial_stock;                   // How much available at the beginning
  int total_avail;                     // How much total in reserve
  int stock;                           // How much do you currently have
  int reserve;                         // How much is still left to system
  int additional_per_day;              // How much can be introduced into the system on a given day
  
  int start_day;                       // When to start production

  int *strains;
  int num_strains;
  
  // for statistics
  int number_delivered;
  int number_effective;
  
protected:
  Vaccine() { }
};

#endif
