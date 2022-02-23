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
// File: Vaccine_Priority_Policies.cc
//

#include <iostream>
#include <vector>

#include "Decision.h"
#include "Vaccine_Priority_Policies.h"
#include "Vaccine_Priority_Decisions.h"
#include "Policy.h"
#include "Manager.h"
#include "Vaccine_Manager.h"

Vaccine_Priority_Policy_No_Priority::Vaccine_Priority_Policy_No_Priority(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  
  Name = "Vaccine Priority Policy - No Priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_No_Priority(this));
}

Vaccine_Priority_Policy_Specific_Age::Vaccine_Priority_Policy_Specific_Age(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  
  Name = "Vaccine Priority Policy - Sepcific Age Group";
  decision_list.push_back(new Vaccine_Priority_Decision_Specific_Age(this));
}

Vaccine_Priority_Policy_ACIP::Vaccine_Priority_Policy_ACIP(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  
  Name = "Vaccine Priority Policy - ACIP Priority";
  decision_list.push_back(new Vaccine_Priority_Decision_Specific_Age(this));
  decision_list.push_back(new Vaccine_Priority_Decision_Pregnant(this));
  decision_list.push_back(new Vaccine_Priority_Decision_At_Risk(this));
}

Vaccine_Priority_Policy_Phase_Age::Vaccine_Priority_Policy_Phase_Age(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - Age specific phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_Age(this));
  
}
  
Vaccine_Priority_Policy_Phase_Essential_Workers::Vaccine_Priority_Policy_Phase_Essential_Workers(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - Essential workers phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_Essential_Workers(this));
  
}

Vaccine_Priority_Policy_Phase_Teachers::Vaccine_Priority_Policy_Phase_Teachers(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - Teachers phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_Teachers(this));
  
}

Vaccine_Priority_Policy_Phase_Comorbidities::Vaccine_Priority_Policy_Phase_Comorbidities(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - Comorbidities phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_Comorbidities(this));
  
}

Vaccine_Priority_Policy_Phase_LTC::Vaccine_Priority_Policy_Phase_LTC(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - Long term care facilities phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_LTC(this));
  
}

Vaccine_Priority_Policy_Phase_No_Priority::Vaccine_Priority_Policy_Phase_No_Priority(Vaccine_Manager *vcm):
  Policy(dynamic_cast <Manager *> (vcm)){
  Name = "Vaccine Priority Policy - No priority in phased priority";
  
  decision_list.push_back(new Vaccine_Priority_Decision_Phase_No_Priority(this));
  
}

