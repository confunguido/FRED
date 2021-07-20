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
// File: Vaccine.cc
//
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "Vaccine.h"
#include "Vaccine_Dose.h"
#include "Global.h"
#include "Disease.h"
#include "Disease_List.h"
Vaccine::Vaccine(string _name, int _id, int _disease, 
                 int _total_avail, int _additional_per_day, 
                 int _start_day, int num_strains, int *_strains){
  name =               _name;
  id =                 _id;
  disease =             _disease;
  additional_per_day = _additional_per_day;
  start_day =          _start_day;
  strains =            _strains;
  
  initial_stock = 0;
  stock = 0;
  reserve = _total_avail;
  total_avail = _total_avail;
  number_delivered = 0;
  number_effective = 0;
  disease_specific_efficacy_modifier.clear();
  disease_specific_efficacy_symp_modifier.clear();
  disease_specific_efficacy_hosp_modifier.clear();
}

Vaccine::~Vaccine(){ 
  for(unsigned int i = 0; i < doses.size(); i++) delete doses[i];
  delete strains;
}

void Vaccine::add_dose(Vaccine_Dose* _vaccine_dose) {
  doses.push_back(_vaccine_dose);
}

void Vaccine::print() const {
  cout << "Name = \t\t\t\t" <<name << "\n";
  cout << "Applied to disease = \t\t" << disease << "\n";
  cout << "Initial Stock = \t\t" << initial_stock << "\n";
  cout << "Total Available = \t\t"<< total_avail << "\n";
  cout << "Amount left to system = \t" << reserve << "\n";
  cout << "Additional Stock per day =\t" << additional_per_day << "\n";
  cout << "Starting on day = \t\t" << start_day << "\n";
  cout << "Dose Information\n";
  for(unsigned int i=0;i<doses.size();i++){
    cout <<"Dose #"<<i+1 << "\n";
    doses[i]->print();
  }
  if(Global::Enable_Disease_Cross_Protection == true && Global::Diseases.get_number_of_diseases() > 1){
    cout << "Differential efficacy information for " << this->disease_specific_efficacy_modifier.size() <<" diseases\n";
    for(int dis_id = 0; dis_id < Global::Diseases.get_number_of_diseases(); ++dis_id){
      cout << "Disease " << dis_id << "Eff. modifier: " << disease_specific_efficacy_modifier[dis_id] << "\n";
      cout << "Disease " << dis_id << "Eff. Symp. modifier: " << disease_specific_efficacy_symp_modifier[dis_id] << "\n";
      cout << "Disease " << dis_id << "Eff. Hosp. modifier: " << disease_specific_efficacy_hosp_modifier[dis_id] << "\n";
    }
  }
}

void Vaccine::reset() {
  stock = 0;
  reserve = total_avail;
}

void Vaccine::update(int day) {
  if(day >= start_day) add_stock(additional_per_day);
}

void Vaccine::update(int day, int add_vac) {
  if(day >= start_day) add_stock(add_vac);
}

int Vaccine::get_strain(int i) {
  if(i < num_strains) return strains[i];
  else return -1;
}
