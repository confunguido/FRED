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
// File: Vaccine_Manager.h
//

#ifndef _FRED_VACCINE_MANAGER_H
#define _FRED_VACCINE_MANAGER_H

#define VACC_NO_PRIORITY  0
#define VACC_AGE_PRIORITY 1
#define VACC_ACIP_PRIORITY 2

#define VACC_DOSE_NO_PRIORITY 0
#define VACC_DOSE_FIRST_PRIORITY 1
#define VACC_DOSE_RAND_PRIORITY 2
#define VACC_DOSE_LAST_PRIORITY 3
#define VACC_DOSE_SEPARATE_PRIORITY 4

#include <list>
#include <vector>
#include <string>
#include "Manager.h"

using namespace std;

class Manager;
class Population;
class Vaccines;
class Person;
class Policy;
class Timestep_Map;
class Events;

class Vaccine_Manager: public Manager {
  //Vaccine_Manager handles a stock of vaccines
public:
  // Construction and Destruction
  Vaccine_Manager();
  Vaccine_Manager(Population * _pop);
  ~Vaccine_Manager();
  
  //Parameters Access
  bool do_vaccination() const {
    return this->do_vacc;
  }

  Vaccines * get_vaccines() const {
    return this->vaccine_package;
  }

  list<Person*> get_priority_queue() const {
    return this->priority_queue;
  }

  list<Person*> get_queue() const {
    return this->queue;}
  int get_number_in_priority_queue() const {
    return this->priority_queue.size();
  }

  int get_number_in_reg_queue() const {
    return this->queue.size();
  }

  int get_current_vaccine_capacity() const {
    return this->current_vaccine_capacity;
  }
  
  // Vaccination Specific Procedures
  void fill_queues();
  void vaccinate(int day);
  void add_to_queue(Person * person);                 //Adds person to queue based on current policies
  void remove_from_queue(Person * person);            //Remove person from the vaccine queue
  void add_to_priority_queue_random(Person * person); //Adds person to the priority queue in a random spot
  void add_to_regular_queue_random(Person * person);  //Adds person to the regular queue in a random spot
  void add_to_priority_queue_begin(Person * person);  //Adds person to the beginning of the priority queue
  void add_to_priority_queue_end(Person * person);    //Adds person to the end of the priority queue
  void add_to_next_dose_queue_end(Person * person); // Adds person to a different priority queue for doses
  void add_next_priority_to_queues();
  
  //Paramters Access Members
  int get_vaccine_priority_age_low() const {
    return vaccine_priority_age_low;
  }

  int get_vaccine_priority_age_high() const {
    return this->vaccine_priority_age_high;
  }

  int get_vaccine_priority_age_low(int n) const {
    if(vaccine_priority_phases_age_low.size() > n){
      return vaccine_priority_phases_age_low[n];
    }else{
      return -120;
    }
  }

  int get_vaccine_priority_age_high(int n) const {
    if(vaccine_priority_phases_age_high.size() > n){
      return this->vaccine_priority_phases_age_high[n];
    }else{
      return -1;
    }       
  }

  double get_vaccine_priority_pop_prob(int n) const {
    if(vaccine_priority_phases_pop_prob.size() > n){
      return this->vaccine_priority_phases_pop_prob[n];
    }else{
      return 0.0;
    }       
  }

  int get_vaccine_dose_priority() const {
    return this->vaccine_dose_priority;
  }

  string get_vaccine_dose_priority_string() const;
  
  // Utility Members
  void update(int day);
  void reset();
  void print();
  
private:
  // events processing
  void process_vaccine_immunity_start_events(int day);
  void process_vaccine_next_dose_events(int day);
  void process_vaccine_immunity_end_events(int day);
  Events* vaccine_next_dose_event_queue;
  Events* vaccine_immunity_start_event_queue;  
  Events* vaccine_immunity_end_event_queue;
  
  Vaccines* vaccine_package;             //Pointer to the vaccines that this manager oversees
  list<Person *> priority_queue;         //Queue for the priority agents
  list<Person *> queue;                  //Queue for everyone else
  list<Person *> next_dose_queue;

  double vaccine_acceptance_prob;
  
  //Parameters from Input 
  bool do_vacc;                           //Is Vaccination being performed
  bool vaccine_priority_only;             //True - Vaccinate only the priority
  bool vaccinate_symptomatics;            //True - Include people that have had symptoms in the vaccination
  bool refresh_vaccine_queues_daily;      //True - people queue up in random order each day
  bool enable_vaccine_priority_discrete_refill;
  
  int vaccine_priority_age_low;           //Age specific priority
  int vaccine_priority_age_high;
  vector<int>vaccine_priority_phases_age_low;
  vector<int>vaccine_priority_phases_age_high;
  vector<int>vaccine_priority_phases_ID;
  vector<int>vaccine_priority_timing_vector;
  vector<double>vaccine_priority_phases_pop_prob;
  std::vector<std::list<Person*>>priority_queue_vector;

  int current_priority_included;
  int vaccine_dose_priority;              //Defines where people getting multiple doses fit in the queue
                                          // See defines above for values
  
  Timestep_Map * vaccination_capacity_map; // How many people can be vaccinated now,
  // gets its value from the capacity change list
  int current_vaccine_capacity;           // variable to keep track of how many persons this 
                                          // can vaccinate each timestep.
};


#endif


