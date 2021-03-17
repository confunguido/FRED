/*
  This file is part of the FRED system.

  Copyright (c) 2010-2015, University of Pittsburgh, John Grefenstette,
  Shawn Brown, Roni Rosenfield, Alona Fyshe, David Galloway, Nathan
  Stone, Jay DePasse, Anuroop Sriram, and Donald Burke.

  Licensed under the BSD 3-Clause license.  See the file "LICENSE" for
  more information.
*/

//
// File: Vaccine_Priority_Decision.h
//

#ifndef _FRED_VACCINE_PRIORITY_DECISIONS_H
#define _FRED_VACCINE_PRIORITY_DECISIONS_H


#include "Decision.h"

class Policy;
class Person;

class Vaccine_Priority_Decision_Specific_Age: public Decision {
public:
  Vaccine_Priority_Decision_Specific_Age(Policy* p);
  Vaccine_Priority_Decision_Specific_Age();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Pregnant: public Decision {
public:
  Vaccine_Priority_Decision_Pregnant(Policy* p);
  Vaccine_Priority_Decision_Pregnant();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_At_Risk: public Decision {
public:
  Vaccine_Priority_Decision_At_Risk(Policy* p);
  Vaccine_Priority_Decision_At_Risk();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_No_Priority: public Decision {
public:
  Vaccine_Priority_Decision_No_Priority(Policy *p);
  Vaccine_Priority_Decision_No_Priority();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_Age: public Decision {
public:
  Vaccine_Priority_Decision_Phase_Age(Policy *p);
  Vaccine_Priority_Decision_Phase_Age();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_Essential_Workers: public Decision {
public:
  Vaccine_Priority_Decision_Phase_Essential_Workers(Policy *p);
  Vaccine_Priority_Decision_Phase_Essential_Workers();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_Teachers: public Decision {
public:
  Vaccine_Priority_Decision_Phase_Teachers(Policy *p);
  Vaccine_Priority_Decision_Phase_Teachers();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_Comorbidities: public Decision {
public:
  Vaccine_Priority_Decision_Phase_Comorbidities(Policy *p);
  Vaccine_Priority_Decision_Phase_Comorbidities();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_LTC: public Decision {
public:
  Vaccine_Priority_Decision_Phase_LTC(Policy *p);
  Vaccine_Priority_Decision_Phase_LTC();
  int evaluate(Person* person, int disease, int day);
};

class Vaccine_Priority_Decision_Phase_No_Priority: public Decision {
public:
  Vaccine_Priority_Decision_Phase_No_Priority(Policy *p);
  Vaccine_Priority_Decision_Phase_No_Priority();
  int evaluate(Person* person, int disease, int day);
};

#endif
