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
// File: Vaccine_Priority_Policies.cpp
//

#include "Vaccine_Priority_Decisions.h"
#include "Vaccine_Priority_Policies.h"
#include "Manager.h"
#include "Vaccine_Manager.h"
#include "Vaccines.h"
#include "Random.h"
#include "Person.h"
#include "Health.h"
#include "Demographics.h"
#include <iostream>
#include <string>

Vaccine_Priority_Decision_Specific_Age::Vaccine_Priority_Decision_Specific_Age() : Decision(){ }

Vaccine_Priority_Decision_Specific_Age::Vaccine_Priority_Decision_Specific_Age(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Specific Age";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Specific_Age::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int low_age = vcm->get_vaccine_priority_age_low();
  int high_age = vcm->get_vaccine_priority_age_high();
  
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    return 1;
  }
  return -1;
}

Vaccine_Priority_Decision_Pregnant::Vaccine_Priority_Decision_Pregnant() : Decision(){ }

Vaccine_Priority_Decision_Pregnant::Vaccine_Priority_Decision_Pregnant(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision to Include Pregnant Women";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Pregnant::evaluate(Person* person, int disease, int day){
  if(person->get_demographics()->is_pregnant()) return 1;
  return -1;
}

Vaccine_Priority_Decision_At_Risk::Vaccine_Priority_Decision_At_Risk() : Decision(){ }

Vaccine_Priority_Decision_At_Risk::Vaccine_Priority_Decision_At_Risk(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision to Include At_Risk";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_At_Risk::evaluate(Person* person, int disease, int day){
  if(person->get_health()->is_at_risk(disease)) return 1;
  return -1;
}


Vaccine_Priority_Decision_No_Priority::Vaccine_Priority_Decision_No_Priority() : Decision() { }

Vaccine_Priority_Decision_No_Priority::Vaccine_Priority_Decision_No_Priority(Policy *p) : Decision(p){
  name = "Vaccine Priority Decision No Priority";
  type = "Y/N";
  policy=p;
}

int Vaccine_Priority_Decision_No_Priority::evaluate(Person* person, int disease, int day){  
  return -1;
}

Vaccine_Priority_Decision_Phase_Age::Vaccine_Priority_Decision_Phase_Age() : Decision(){ }

Vaccine_Priority_Decision_Phase_Age::Vaccine_Priority_Decision_Phase_Age(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases Age";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_Age::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int num_priority = vcm->get_current_policy();
  int low_age = vcm->get_vaccine_priority_age_low(num_priority);
  int high_age = vcm->get_vaccine_priority_age_high(num_priority);
  
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    return 1;
  }
  return -1;
}

Vaccine_Priority_Decision_Phase_Essential_Workers::Vaccine_Priority_Decision_Phase_Essential_Workers() : Decision(){ }

Vaccine_Priority_Decision_Phase_Essential_Workers::Vaccine_Priority_Decision_Phase_Essential_Workers(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases Essential Workers";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_Essential_Workers::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int num_priority = vcm->get_current_policy();
  int low_age = vcm->get_vaccine_priority_age_low(num_priority);
  int high_age = vcm->get_vaccine_priority_age_high(num_priority);
  double pop_prob = vcm->get_vaccine_priority_pop_prob(num_priority);
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    if(Random::draw_random() < pop_prob){
      return 1;
    }else{
      return -1;
    }
  }
  return -1;
}


Vaccine_Priority_Decision_Phase_Teachers::Vaccine_Priority_Decision_Phase_Teachers() : Decision(){ }

Vaccine_Priority_Decision_Phase_Teachers::Vaccine_Priority_Decision_Phase_Teachers(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases Teachers";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_Teachers::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int num_priority = vcm->get_current_policy();
  int low_age = vcm->get_vaccine_priority_age_low(num_priority);
  int high_age = vcm->get_vaccine_priority_age_high(num_priority);
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    if(person->is_teacher()){
      return 1;
    }else{
      return -1;
    }
  }
  return -1;
}

Vaccine_Priority_Decision_Phase_Comorbidities::Vaccine_Priority_Decision_Phase_Comorbidities() : Decision(){ }

Vaccine_Priority_Decision_Phase_Comorbidities::Vaccine_Priority_Decision_Phase_Comorbidities(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases Comorbidities";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_Comorbidities::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int num_priority = vcm->get_current_policy();
  int low_age = vcm->get_vaccine_priority_age_low(num_priority);
  int high_age = vcm->get_vaccine_priority_age_high(num_priority);
  
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    if(person->has_chronic_condition()){
      return 1;
    }else{
      return -1;
    }
  }
  return -1;
}

Vaccine_Priority_Decision_Phase_LTC::Vaccine_Priority_Decision_Phase_LTC() : Decision(){ }

Vaccine_Priority_Decision_Phase_LTC::Vaccine_Priority_Decision_Phase_LTC(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases LTC";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_LTC::evaluate(Person* person, int disease, int day){
  Vaccine_Manager* vcm = dynamic_cast < Vaccine_Manager* >(policy->get_manager());
  int num_priority = vcm->get_current_policy();
  int low_age = vcm->get_vaccine_priority_age_low(num_priority);
  int high_age = vcm->get_vaccine_priority_age_high(num_priority);
  
  if(person->get_age() >= low_age && person->get_age() <= high_age){
    if(person->is_nursing_home_resident()){
      return 1;
    }else{
      return -1;
    }
  }
  return -1;
}

Vaccine_Priority_Decision_Phase_No_Priority::Vaccine_Priority_Decision_Phase_No_Priority() : Decision(){ }

Vaccine_Priority_Decision_Phase_No_Priority::Vaccine_Priority_Decision_Phase_No_Priority(Policy *p): Decision(p) {
  name = "Vaccine Priority Decision Phases No Priority";
  type = "Y/N";
  policy = p;
}

int Vaccine_Priority_Decision_Phase_No_Priority::evaluate(Person* person, int disease, int day){
  return -1;
}
