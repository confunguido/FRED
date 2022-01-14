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
// File: Vaccines.h
//
#ifndef _FRED_VACCINES_H
#define _FRED_VACCINES_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

class Vaccine;
class Vaccine_Dose;

struct Time_Step_Map_Vaccine_Additional {
  int sim_day_start;
  int sim_day_end;
  int additional_vaccine;
  int vaccine_id;
  const std::string to_string() const {
    std::stringstream ss;
    ss << "Time Step Map ";
    ss << " sim_day_start " << sim_day_start;
    ss << " sim_day_end " << sim_day_end;
    ss << " additional stock " << additional_vaccine;
    ss << " vaccine id " << vaccine_id;
    ss << std::endl;
    return ss.str();
  }
};

class Vaccines {
  // Vaccines is a class used to describe a group of Vaccine Classes
public:
  // Creation Operations
  Vaccines() { }
  void setup();
  
  Vaccine *get_vaccine(int i) const { return vaccines[i];}
  
  vector <int> which_vaccines_applicable(double real_age) const;
  int pick_from_applicable_vaccines(double real_age) const;
  int get_total_vaccines_avail_today() const;
  unsigned int get_number_of_vaccine_types() const {return vaccines.size();}
  
  //utility Functions
  void print() const;
  void print_current_stocks() const;
  void update(int day);
  void reset();
private:
  vector < Vaccine* > vaccines;
  static std::vector<Time_Step_Map_Vaccine_Additional*> vaccine_additional_timestep;
}; 

#endif
