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
// File: VaccineDose.cc
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include "Vaccine_Dose.h"
#include "Random.h"

Vaccine_Dose::Vaccine_Dose(Age_Map* _efficacy, Age_Map* _efficacy_symp,
			   Age_Map* _efficacy_hosp, Age_Map* _efficacy_delay,
			   Age_Map* _efficacy_duration, int _days_between_doses,
			   int _mix_match){
  efficacy = _efficacy;
  efficacy_symp = _efficacy_symp;
  efficacy_hosp = _efficacy_hosp;
  efficacy_delay = _efficacy_delay;
  efficacy_duration = _efficacy_duration;
  days_between_doses = _days_between_doses;
  mix_and_match_next_dose = _mix_match;
}

Vaccine_Dose::~Vaccine_Dose(){
  if (efficacy) delete efficacy;
  if (efficacy_symp) delete efficacy_symp;
  if (efficacy_hosp) delete efficacy_hosp;
  if (efficacy_delay) delete efficacy_delay;
}

void Vaccine_Dose::print() const {
  cout << "Time Between Doses:\t " << days_between_doses << "\n";
  cout << "Mix and match next dose:\t " << mix_and_match_next_dose << "\n";
  efficacy->print();
  efficacy_symp->print();
  efficacy_hosp->print();
  efficacy_delay->print();
  efficacy_duration->print();
}

bool Vaccine_Dose::is_within_age(double real_age) const {
  // TODO: This scenario should also consider efficacy against symptoms
  double eff = efficacy->find_value(real_age);
  // printf("age = %.1f  eff = %f\n", real_age, eff);
  if(eff != 0.0){
    return true;
  }
  return false;
}


double Vaccine_Dose::get_duration_of_immunity(double real_age) {
  double expected_duration = efficacy_duration->find_value(real_age);
  // select a value from an exponential distribution with mean expected_duration
  double actual_duration = 0.0;
  if (expected_duration > 0.0) {
    actual_duration = Random::draw_exponential(1.0 / expected_duration);
  }
  return actual_duration;
}
