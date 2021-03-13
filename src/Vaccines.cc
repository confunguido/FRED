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
// File: Vaccines.cc
//


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

#include "Vaccines.h"
#include "Vaccine.h"
#include "Vaccine_Dose.h"
#include "Random.h"
#include "Age_Map.h"

std::vector<Time_Step_Map_Vaccine_Additional*> Vaccines::vaccine_additional_timestep;

void Vaccines::setup(void) {

  int number_vacc;
  Params::get_param_from_string("number_of_vaccines",&number_vacc);
  
  for(int iv=0;iv<number_vacc;iv++) {
    int ta;
    int apd;
    int std;
    int tbd;
    int num_doses;
   
    Params::get_indexed_param("vaccine_number_of_doses",iv,&num_doses);
    Params::get_indexed_param("vaccine_total_avail",iv,&ta);
    Params::get_indexed_param("vaccine_additional_per_day",iv,&apd);
    Params::get_indexed_param("vaccine_starting_day",iv,&std);
    Age_Map* efficacy_duration_map = new Age_Map("Vaccine Efficacy Duration");
    efficacy_duration_map->read_from_input("vaccine_efficacy_duration",iv);

    int nstrains;
    Params::get_indexed_param((char *)"vaccine_strains", iv, &nstrains);
    int *strains = new int[nstrains];
    Params::get_indexed_param_vector<int>("vaccine_strains", iv, strains);

    stringstream name;
    name << "Vaccine#"<<iv+1;
    vaccines.push_back(new Vaccine(name.str(),iv,0,ta,apd,std,nstrains,strains));
    
    for(int id=0;id<num_doses;id++) {
      Age_Map* efficacy_map = new Age_Map("Dose Efficacy");
      Age_Map* efficacy_symp_map = new Age_Map("Dose Efficacy Symptoms");
      Age_Map* efficacy_hosp_map = new Age_Map("Dose Efficacy Hospitalization");      
      Age_Map* efficacy_delay_map = new Age_Map("Dose Efficacy Delay");
      Params::get_double_indexed_param("vaccine_next_dosage_day",iv,id,&tbd);
      efficacy_map->read_from_input("vaccine_dose_efficacy",iv,id);
      efficacy_symp_map->read_from_input("vaccine_dose_efficacy_symptoms",iv,id);
      efficacy_hosp_map->read_from_input("vaccine_dose_efficacy_hospitalization",iv,id);
      efficacy_delay_map->read_from_input("vaccine_dose_efficacy_delay",iv,id);
      vaccines[iv]->add_dose(new Vaccine_Dose(efficacy_map,efficacy_symp_map,efficacy_hosp_map,efficacy_delay_map,efficacy_duration_map,tbd));
    }
  }
  
  if(Global::Enable_Vaccine_Stock_Timeseries_File == true) {
    char map_file_name[FRED_STRING_SIZE];
    Params::get_param_from_string("vaccine_stock_timeseries_file", map_file_name);
    printf("READING Vaccine stock timeseries %s\n", map_file_name);

    // If this parameter is "none", then there is no map
    if(strncmp(map_file_name, "none", 4) != 0){
      Utils::get_fred_file_name(map_file_name);
      printf("READING: vaccine stock timeseries file: %s\n", map_file_name);
      ifstream* ts_input = new ifstream(map_file_name);
      if(!ts_input->is_open()) {
	Utils::fred_abort("Help!  Can't read %s Timestep Map\n", map_file_name);
	abort();
      }
      string line;
      while(getline(*ts_input,line)){
	if(line[0] == '\n' || line[0] == '#') { // empty line or comment
	  continue;
	}
	char cstr[FRED_STRING_SIZE];

	std::strcpy(cstr, line.c_str());
	Time_Step_Map_Vaccine_Additional * tmap = new Time_Step_Map_Vaccine_Additional;
	int n = sscanf(cstr,
		       "%d %d %d %d",
		       &tmap->sim_day_start, &tmap->sim_day_end, &tmap->additional_vaccine, &tmap->vaccine_id);
	printf("LINES: %d\n",n);
	if(n < 4) {
	  Utils::fred_abort("Need to specify at least SimulationDayStart, SimulationDayEnd, AdditionalStock, and vaccine ID.");
	}
	if(tmap->vaccine_id >= 0 && tmap->additional_vaccine >= 0){
	  this->vaccine_additional_timestep.push_back(tmap);
	}
      }
      ts_input->close();
    }
    
    for(int i = 0; i < this->vaccine_additional_timestep.size(); ++i){
      string ss = this->vaccine_additional_timestep[i]->to_string();
      printf("%s\n", ss.c_str());
    }
    
    printf("Finished with reading vaccine stock timeseries file\n");
  }
}  

void Vaccines::print() const {
  cout <<"Vaccine Package Information\n";
  cout <<"There are "<<vaccines.size() <<" vaccines in the package\n";
  fflush(stdout);
  for(unsigned int i=0;i<vaccines.size(); i++){
    vaccines[i]->print();
  }
  fflush(stdout);
}

void Vaccines::print_current_stocks() const {
  cout << "Vaccine Stockk Information\n";
  cout << "\nVaccines# " << "Current Stock      " << "Current Reserve    \n";
  for(unsigned int i=0; i<vaccines.size(); i++) {
    cout << setw(10) << i+1 << setw(20) << vaccines[i]->get_current_stock()
	 << setw(20) << vaccines[i]->get_current_reserve() << "\n";
  }
}

void Vaccines::update(int day) {
  for(unsigned int i=0;i<vaccines.size(); i++) {
    if(Global::Enable_Vaccine_Stock_Timeseries_File == false){
      vaccines[i]->update(day);
    }else{
      // find vaccine ID and day
      int additional_stock_today = 0;
      for(int j = 0; j < this->vaccine_additional_timestep.size(); ++j){
	Time_Step_Map_Vaccine_Additional* tmap = this->vaccine_additional_timestep[j];
	if(tmap->sim_day_start <= day && day <= tmap->sim_day_end && tmap->vaccine_id == i) {
	  additional_stock_today = tmap->additional_vaccine;	  
	}
      }
      vaccines[i]->update(day, additional_stock_today);
    }
  }
}

void Vaccines::reset() {
  for(unsigned int i=0;i<vaccines.size();i++) {
    vaccines[i]->reset();
  }
}


int Vaccines::pick_from_applicable_vaccines(double real_age) const {
  vector <int> app_vaccs;
  for(unsigned int i=0;i<vaccines.size();i++){
    // if first dose is applicable, add to vector.
    if(vaccines[i]->get_dose(0)->is_within_age(real_age) &&
       vaccines[i]->get_current_stock() > 0){
      app_vaccs.push_back(i);
    }
  }
  
  if(app_vaccs.size() == 0){ return -1; }
  
  int randnum = 0;
  if(app_vaccs.size() > 1){
    randnum = (int)(Random::draw_random()*app_vaccs.size());
  }
  return app_vaccs[randnum];
}

int Vaccines::get_total_vaccines_avail_today() const {
  int total=0;
  for(unsigned int i=0;i<vaccines.size();i++){
    total += vaccines[i]->get_current_stock();
  }
  return total;
}
