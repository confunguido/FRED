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
// File: Vaccine_Manager.cpp
//
#include "Manager.h"
#include "Vaccine_Manager.h"
#include "Policy.h"
#include "Vaccine_Priority_Policies.h"
#include "Population.h"
#include "Vaccines.h"
#include "Vaccine.h"
#include "Person.h"
#include "Health.h"
#include "Activities.h"
#include "Params.h"
#include "Random.h"
#include "Global.h"
#include "Timestep_Map.h"
#include "Utils.h"
#include "Tracker.h"
#include "Events.h"

#include <algorithm>

Vaccine_Manager::Vaccine_Manager() {
  this->vaccine_package = NULL;
  this->vaccine_priority_age_low = -1;
  this->vaccine_priority_age_high = -1;
  this->vaccine_priority_timing_vector.clear();
  this->vaccine_priority_phases_age_low.clear();
  this->vaccine_priority_phases_age_high.clear();
  this->vaccine_priority_phases_ID.clear();
  this->vaccine_priority_phases_pop_prob.clear();
  this->priority_queue_vector.clear();
  this->current_vaccine_capacity = -1;
  this->current_priority_included = -1;
  this->vaccine_priority_only = false;
  this->enable_vaccine_priority_discrete_refill = false;
  this->vaccination_capacity_map = NULL;
  this->vaccine_acceptance_prob = 1.0;
  this->do_vacc = false;
  this->vaccine_dose_priority = -1;
  this->refresh_vaccine_queues_daily = false;
  this->vaccinate_symptomatics = false;
  this->vaccine_next_dose_event_queue = new Events();
  this->vaccine_immunity_start_event_queue = new Events();
  this->vaccine_immunity_end_event_queue = new Events();
  printf("Finished creating vaccine manager ()\n");
}

Vaccine_Manager::Vaccine_Manager(Population *_pop) :
  Manager(_pop) {
  printf("Vaccine manager entered\n");
  this->priority_queue_vector.clear();
  this->pop = _pop;

  this->vaccine_next_dose_event_queue = new Events();
  this->vaccine_immunity_start_event_queue = new Events();
  this->vaccine_immunity_end_event_queue = new Events();
  printf("Created queues in vaccine manager (pop)\n");

  this->vaccine_package = new Vaccines();
  int num_vaccs = 0;
  Params::get_param_from_string("number_of_vaccines", &num_vaccs);
  if(num_vaccs > 0) {
    printf("Print vaccines from the package\n");
    this->vaccine_package->setup();
    this->vaccine_package->print();
    this->do_vacc = true;
  } else {    // No vaccination specified.
    this->vaccine_priority_age_low = -1;
    this->vaccine_priority_age_high = -1;
    this->vaccination_capacity_map = NULL;
    this->current_vaccine_capacity = -1;
    this->vaccine_dose_priority = -1;
    this->vaccine_priority_only = false;
    this->enable_vaccine_priority_discrete_refill = false;
    this->current_priority_included = -1;
    this->do_vacc = false;
    this->priority_queue_vector.clear();
    return;
  }

  if(Global::Enable_Vaccination_Phases == false){
    // ACIP Priority takes precidence
    int do_acip_priority;
    current_policy = VACC_NO_PRIORITY;
    Params::get_param_from_string("vaccine_prioritize_acip", &do_acip_priority);
    if(do_acip_priority == 1) {
      cout << "Vaccination Priority using ACIP recommendations\n";
      cout << "   Includes: \n";
      cout << "        Ages 0 to 24\n";
      cout << "        Pregnant Women\n";
      cout << "        Persons at risk for complications\n";
      this->current_policy = VACC_ACIP_PRIORITY;
      this->vaccine_priority_age_low = 0;
      this->vaccine_priority_age_high = 24;
    } else {
      int do_age_priority;
      Params::get_param_from_string("vaccine_prioritize_by_age", &do_age_priority);
      if(do_age_priority) {
	cout << "Vaccination Priority by Age\n";
	this->current_policy = VACC_AGE_PRIORITY;
	Params::get_param_from_string("vaccine_priority_age_low", &this->vaccine_priority_age_low);
	Params::get_param_from_string("vaccine_priority_age_high", &this->vaccine_priority_age_high);
	cout << "      Between Ages " << this->vaccine_priority_age_low << " and "
	     << this->vaccine_priority_age_high << "\n";
      } else {
	this->vaccine_priority_age_low = 0;
	this->vaccine_priority_age_high = Demographics::MAX_AGE;
      }
    }
  }

  // Vaccine acceptance if individual behaviors are disabled
  Params::get_param_from_string("vaccine_acceptance_probability", &this->vaccine_acceptance_prob);
  
  // should we vaccinate anyone outside of the priority class
  int vacc_pri_only;
  this->vaccine_priority_only = false;
  Params::get_param_from_string("vaccine_priority_only", &vacc_pri_only);
  if(vacc_pri_only) {
    this->vaccine_priority_only = true;
    cout << "      Vaccinating only the priority groups\n";
  }

  // should we exclude people that have had symptomatic infections?
  int vacc_sympt_exclude;
  this->vaccinate_symptomatics = false;
  Params::get_param_from_string("vaccinate_symptomatics", &vacc_sympt_exclude);
  if(vacc_sympt_exclude) {
    this->vaccinate_symptomatics = true;
    cout << "      Vaccinating symptomatics\n";
  }

  // get vaccine_dose_priority
  Params::get_param_from_string("vaccine_dose_priority", &this->vaccine_dose_priority);
  assert(this->vaccine_dose_priority < 5);
  //get_param((char*)"vaccination_capacity",&vaccination_capacity);
  this->vaccination_capacity_map = new Timestep_Map("vaccination_capacity");
  this->vaccination_capacity_map->read_map();
  if(Global::Verbose > 1) {
    this->vaccination_capacity_map->print();
  }

  int refresh;
  Params::get_param_from_string("refresh_vaccine_queues_daily", &refresh);
  this->refresh_vaccine_queues_daily = (refresh > 0);

  if(Global::Enable_Vaccination_Phases == true){
    // Fill the vaccine manager policies based on the number of phases and each phase identifier
    int number_of_phases, N_phases_age_low, N_phases_age_high, N_phases_id,N_phases_pop_prob,N_phase_timing,tmp_int;
    Params::get_param_from_string("vaccination_phases_enable_discrete_timing", &(tmp_int));
    this->enable_vaccine_priority_discrete_refill = (tmp_int == 1 ? true : false);
    if(this->enable_vaccine_priority_discrete_refill == true){
      Params::get_param_from_string("vaccination_phases_discrete_timing", &(N_phase_timing));
    }
    Params::get_param_from_string("vaccination_phases_names", &(number_of_phases));
    Params::get_param_from_string("vaccination_phases_age_low", &(N_phases_age_low));
    Params::get_param_from_string("vaccination_phases_age_high", &(N_phases_age_high));
    Params::get_param_from_string("vaccination_phases_id", &(N_phases_id));
    Params::get_param_from_string("vaccination_phases_pop_prob", &(N_phases_pop_prob));
    printf("Vaccination phases enabled:: N_names:%d, N_age_low: %d, N_age_high: %d, N_phases_id %d, N_pop_prob %d\n", number_of_phases,N_phases_age_low, N_phases_age_high, N_phases_id, N_phases_pop_prob);
    if(!(number_of_phases == N_phases_age_low && number_of_phases == N_phases_age_high && number_of_phases == N_phases_id && number_of_phases == N_phases_pop_prob)){
      Utils::fred_abort("Number of phases should be the same as the ages low and high for vaccination phases, and their IDs, as well as their pop_prob even if pop_prob is 0");
    }
    if(this->enable_vaccine_priority_discrete_refill == true){
      if(N_phase_timing != number_of_phases){
	this->enable_vaccine_priority_discrete_refill = false;
      }
    }
    // If phases not specified, create a non-priority one
    if(number_of_phases <= 0){
      // Make sure this would work
      current_policy = VACC_NO_PRIORITY;
      this->policies.push_back(new Vaccine_Priority_Policy_No_Priority(this));
      this->vaccine_priority_age_low = 0;
      this->vaccine_priority_age_high = Demographics::MAX_AGE;
      printf("Vaccination phases enabled, but phases not specified. Assuming no priority\n");
    }else{
      std::vector<string> vaccination_phase_names_un;
      std::vector<string> vaccination_phase_names;
      std::vector<double>vaccine_priority_phases_age_low_un;
      std::vector<double>vaccine_priority_phases_timing_un;
      std::vector<double>vaccine_priority_phases_age_high_un;
      std::vector<double>vaccine_priority_phases_ID_un;
      std::vector<double>vaccine_priority_phases_pop_prob_un;
      Params::get_param_vector_from_string((char*)"vaccination_phases_names", vaccination_phase_names_un);
      char age_low_str[MAX_PARAM_SIZE], age_high_str[MAX_PARAM_SIZE], phase_id_str[MAX_PARAM_SIZE], phase_pop_str[MAX_PARAM_SIZE], phase_time_str[MAX_PARAM_SIZE];
      Params::get_param((char*)"vaccination_phases_age_low", age_low_str);
      Params::get_param((char*)"vaccination_phases_age_high",age_high_str);
      Params::get_param((char*)"vaccination_phases_id",phase_id_str);
      Params::get_param((char*)"vaccination_phases_pop_prob",phase_pop_str);
      Params::get_param_vector_from_string(age_low_str, vaccine_priority_phases_age_low_un);
      Params::get_param_vector_from_string(age_high_str, vaccine_priority_phases_age_high_un);
      Params::get_param_vector_from_string(phase_id_str, vaccine_priority_phases_ID_un);
      Params::get_param_vector_from_string(phase_pop_str, vaccine_priority_phases_pop_prob_un);
      if(this->enable_vaccine_priority_discrete_refill == true){
	Params::get_param((char*)"vaccination_phases_discrete_timing",phase_time_str);
	Params::get_param_vector_from_string(phase_time_str, vaccine_priority_phases_timing_un);
      }

      std::vector<int>sorted_phase_id;
      std::vector<int>sorted_phase_indices;
      std::set<int>unique_phase_id;
      printf("Trying to reassign vaccine priority phases\n");
      for(int i = 0; i < vaccine_priority_phases_ID_un.size(); ++i){
	sorted_phase_id.push_back((int)vaccine_priority_phases_ID_un[i]);
      }
      printf("Reassigned to sorted phase id\n");
      std::sort(sorted_phase_id.begin(), sorted_phase_id.end());
      for(int i = 0; i < sorted_phase_id.size(); ++i){
	unique_phase_id.insert(sorted_phase_id[i]);
      }
      printf("Assigned phases to a set\n");
      for(auto it : unique_phase_id){
	for(int i = 0; i < vaccine_priority_phases_ID_un.size(); ++i){
	  if(it ==  (int)vaccine_priority_phases_ID_un[i]){
	    sorted_phase_indices.push_back(i);
	    this->vaccine_priority_phases_age_low.push_back((int) vaccine_priority_phases_age_low_un[i]);
	    this->vaccine_priority_phases_age_high.push_back((int)vaccine_priority_phases_age_high_un[i]);
	    this->vaccine_priority_phases_ID.push_back((int) vaccine_priority_phases_ID_un[i]);
	    this->vaccine_priority_phases_pop_prob.push_back(vaccine_priority_phases_pop_prob_un[i]);
	    this->vaccine_priority_timing_vector.push_back(vaccine_priority_phases_timing_un[i]);
	    vaccination_phase_names.push_back(vaccination_phase_names_un[i]);
	  }
	}
      }            
      // Should we first sort and then add to the list of policies?
      for(int pp = 0; pp < this->vaccine_priority_phases_ID.size(); pp++){	
	// Read each priority group and find an appropriate policy, if non existent, add no priority
	printf("index = %d, priority id: %d, priority name %s\n", pp,(int)this->vaccine_priority_phases_ID[pp], vaccination_phase_names[pp].c_str());
	if(vaccination_phase_names[pp] == "age"){
	  printf("Vaccination priority %d, priority based on age. Min Age %d, Max Age %d, ID %d\n", pp, (int)this->vaccine_priority_phases_age_low[pp], (int)this->vaccine_priority_phases_age_high[pp], (int) this->vaccine_priority_phases_ID[pp]);
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_Age(this));
	}else if(vaccination_phase_names[pp] == "essentialworkers"){
	  printf("Vaccination priority %d, priority to essential workers, ID %d\n", pp, (int)this->vaccine_priority_phases_ID[pp]);
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_Essential_Workers(this));
	}else if(vaccination_phase_names[pp] == "comorbidity"){
	  printf("Vaccination priority %d, priority to comorbidity groups ID %d\n", pp,(int)this->vaccine_priority_phases_ID[pp]);
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_Comorbidities(this));
	}else if(vaccination_phase_names[pp] == "ltc"){
	  printf("Vaccination priority %d, priority to LTC groups ID %d\n", pp,(int)this->vaccine_priority_phases_ID[pp]);
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_LTC(this));
	}else if(vaccination_phase_names[pp] == "teachers"){
	  printf("Vaccination priority %d, priority to teachers groups ID %d\n", pp,(int)this->vaccine_priority_phases_ID[pp]);
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_Teachers(this));
	}else{
	  printf("Vaccination priority not available, adding a non-priority phase\n");
	  this->vaccine_priority_phases_age_low[pp] = 0;
	  this->vaccine_priority_phases_age_high[pp] = Demographics::MAX_AGE;
	  this->policies.push_back(new Vaccine_Priority_Policy_Phase_No_Priority(this));
	}
      }
    }
  }else{
    // Need to fill the Vaccine_Manager Policies
    this->policies.push_back(new Vaccine_Priority_Policy_No_Priority(this));
    this->policies.push_back(new Vaccine_Priority_Policy_Specific_Age(this));
    this->policies.push_back(new Vaccine_Priority_Policy_ACIP(this));
  }

  if(Global::Enable_Vaccination_Phases == true && this->vaccine_priority_phases_ID.size() > 0){
    for(int i = 0; i < this->policies.size(); i++){
      this->policies[i]->print();
      printf("POLICY %d ID %d---------\nAge min: %d, Age max: %d, Prop. Pop: %.3f Timing %d\n", i, this->vaccine_priority_phases_ID[i], this->vaccine_priority_phases_age_low[i], this->vaccine_priority_phases_age_high[i],
	     this->vaccine_priority_phases_pop_prob[i],this->vaccine_priority_timing_vector[i]);
    }
  }
};


Vaccine_Manager::~Vaccine_Manager() {
  if(this->vaccine_package != NULL) {
    delete this->vaccine_package;
  }
  if(this->vaccination_capacity_map != NULL) {
    delete this->vaccination_capacity_map;
  }
}

void Vaccine_Manager::add_next_priority_to_queues(){
  if(!this->do_vacc || this->enable_vaccine_priority_discrete_refill == false || Global::Enable_Vaccination_Phases == false || this->vaccine_priority_phases_ID.size() == 0) {
    return;
  }
  
  if(this->current_priority_included >= (this->policies.size() - 1) || this->current_priority_included == -1){
    return;
  }

  // Load the first queue with policy 0
  printf("Re-Filling up vaccination priority queues\n");    
  printf("Priority queue has %lu size, regular queue %lu size\n", this->priority_queue.size(), this->queue.size());
  
  int ii = this->current_priority_included + 1;
  printf("Priority queue for phase %d - %d has size: %lu, regular queue has size %lu\n", (int)this->vaccine_priority_phases_ID[ii], ii, this->priority_queue_vector[ii].size(), this->queue.size());
  
  vector<Person *> rand_priority_queue(this->priority_queue_vector[ii].size());
  copy(this->priority_queue_vector[ii].begin(), this->priority_queue_vector[ii].end(), rand_priority_queue.begin());
  FYShuffle<Person *>(rand_priority_queue);
  copy(rand_priority_queue.begin(), rand_priority_queue.end(),
       this->priority_queue_vector[ii].begin());
      
  while(!this->priority_queue_vector[ii].empty()){
    this->priority_queue.push_back(this->priority_queue_vector[ii].front());
    this->priority_queue_vector[ii].pop_front();
  }

  printf("Priority queue has %lu size, regular queue %lu size\n", this->priority_queue.size(), this->queue.size());    

  vector<Person *> random_queue(this->queue.size());
  copy(this->queue.begin(), this->queue.end(), random_queue.begin());
  FYShuffle<Person *>(random_queue);
  copy(random_queue.begin(), random_queue.end(), this->queue.begin());


  vector<Person *> random_priority_queue(this->priority_queue.size());
  copy(this->priority_queue.begin(), this->priority_queue.end(), random_priority_queue.begin());
  FYShuffle<Person *>(random_priority_queue);
  copy(random_priority_queue.begin(), random_priority_queue.end(),
       this->priority_queue.begin());
  this->current_priority_included++;
  if(Global::Verbose > 0) {
    cout << "After refilling!! Vaccine Queue Stats \n";
    cout << "   Number in Priority Queue      = " << this->priority_queue.size() << "\n";
    cout << "   Number in Regular Queue       = " << this->queue.size() << "\n";
    cout << "   Total Agents in Vaccine Queue = "
         << this->priority_queue.size() + this->queue.size() << "\n";
  }
}
void Vaccine_Manager::fill_queues() {

  if(!this->do_vacc) {
    return;
  }
  /*
    We need to loop over the entire population that the Manager oversees to put them in a queue.
    If phases are enabled, loop over all the policies, in order and then add to the queue
    Maybe check if policies are overlapping, in that case, only enroll in one queue
   */
  if(Global::Enable_Vaccination_Phases == true && this->vaccine_priority_phases_ID.size() > 0){
    // Load the first queue with policy 0 <---- MOVE THIS TO INITIALIZATION
    printf("Filling up vaccination priority queues\n");
    //std::vector<std::list<Person*>>priority_queue_vector;
    
    for(int ii = 0; ii < this->vaccine_priority_phases_ID.size(); ii++){
      std::list<Person*> tmp_list;
      tmp_list.clear();
      this->priority_queue_vector.push_back(tmp_list);
    }
    printf("Priority queue vector has %d elements, and policies %d\n", this->priority_queue_vector.size(), this->policies.size());
      
    priority_queue.clear();
    for(int ip = 0; ip < pop->get_index_size(); ip++) {
      Person * current_person = this->pop->get_person_by_index(ip);
      if (current_person != NULL) {
	bool current_person_priority = false;
	for(int cc = 0; cc < this->policies.size();cc++){
	  this->current_policy = cc;
	  if(this->policies[current_policy]->choose_first_positive(current_person, 0, 0) == true) {
	    this->priority_queue_vector[current_policy].push_back(current_person);
	    current_person_priority = true;
	    break;
	  }
	}
	if(this->vaccine_priority_only == false && current_person_priority == false){
	  this->queue.push_back(current_person);
	}
      }
    }
    printf("Priority queue has %lu size, regular queue %lu size\n", this->priority_queue.size(), this->queue.size());
    int max_phases = this->priority_queue_vector.size();
    if(this->enable_vaccine_priority_discrete_refill == true){
      max_phases = 1;
      this->current_priority_included = 0;
    }
    for(int ii = 0; ii < max_phases; ii++){
      printf("Priority queue for phase %d - %d has size: %lu, regular queue has size %lu\n", (int)this->vaccine_priority_phases_ID[ii], ii, this->priority_queue_vector[ii].size(), this->queue.size());

      vector<Person *> rand_priority_queue(this->priority_queue_vector[ii].size());
      copy(this->priority_queue_vector[ii].begin(), this->priority_queue_vector[ii].end(), rand_priority_queue.begin());
      FYShuffle<Person *>(rand_priority_queue);
      copy(rand_priority_queue.begin(), rand_priority_queue.end(),
	   this->priority_queue_vector[ii].begin());
      
      while(!this->priority_queue_vector[ii].empty()){
	this->priority_queue.push_back(this->priority_queue_vector[ii].front());
	this->priority_queue_vector[ii].pop_front();
      }
    }    
    printf("Priority queue has %lu size, regular queue %lu size\n", this->priority_queue.size(), this->queue.size());    
  }else{
    for(int ip = 0; ip < pop->get_index_size(); ip++) {
      Person * current_person = this->pop->get_person_by_index(ip);
      if (current_person != NULL) {
	if(this->policies[current_policy]->choose_first_positive(current_person, 0, 0) == true) {
	  priority_queue.push_back(current_person);
	} else {
	  if(this->vaccine_priority_only == false)
	    this->queue.push_back(current_person);
	}
      }
    }
  }

  vector<Person *> random_queue(this->queue.size());
  copy(this->queue.begin(), this->queue.end(), random_queue.begin());
  FYShuffle<Person *>(random_queue);
  copy(random_queue.begin(), random_queue.end(), this->queue.begin());

  // Shouldn't shuffle if priority policies are enabled
  if(!(Global::Enable_Vaccination_Phases == true && this->vaccine_priority_phases_ID.size() > 0)){
    vector<Person *> random_priority_queue(this->priority_queue.size());
    copy(this->priority_queue.begin(), this->priority_queue.end(), random_priority_queue.begin());
    FYShuffle<Person *>(random_priority_queue);
    copy(random_priority_queue.begin(), random_priority_queue.end(),
	 this->priority_queue.begin());
  }

  if(Global::Verbose > 0) {
    cout << "Vaccine Queue Stats \n";
    cout << "   Number in Priority Queue      = " << this->priority_queue.size() << "\n";
    cout << "   Number in Regular Queue       = " << this->queue.size() << "\n";
    cout << "   Total Agents in Vaccine Queue = "
         << this->priority_queue.size() + this->queue.size() << "\n";
  }
  next_dose_queue.clear();
  printf("Size of policies %lu size of phases id %lu\n", this->policies.size(), this->priority_queue_vector.size());
}

void Vaccine_Manager::add_to_queue(Person* person) {
  if(this->policies[this->current_policy]->choose_first_positive(person, 0, 0) == true) {
    add_to_priority_queue_random(person);
  } else {
    if(this->vaccine_priority_only == false) {
      add_to_regular_queue_random(person);
    }
  }
}

void Vaccine_Manager::remove_from_queue(Person* person) {
  // remove the person from the queue if they are in there
  list<Person *>::iterator pq = find(this->priority_queue.begin(), this->priority_queue.end(),
				     person);
  if(pq != this->priority_queue.end()) {
    this->priority_queue.erase(pq);
    return;
  }
  pq = find(this->queue.begin(), this->queue.end(), person);
  if(pq != this->queue.end()) {
    this->queue.erase(pq);
  }
}

void Vaccine_Manager::add_to_priority_queue_random(Person* person) {
  // Find a position to put the person in
  int size = this->priority_queue.size();
  int position = (int)(Random::draw_random()*size);

  list<Person*>::iterator pq = this->priority_queue.begin();
  for(int i = 0; i < position; ++i) {
    ++pq;
  }
  this->priority_queue.insert(pq, person);
}

void Vaccine_Manager::add_to_regular_queue_random(Person* person) {
  // Find a position to put the person in
  int size = this->queue.size();
  int position = (int)(Random::draw_random() * size);

  list<Person*>::iterator pq = this->queue.begin();
  for(int i = 0; i < position; ++i) {
    ++pq;
  }
  this->queue.insert(pq, person);
}

void Vaccine_Manager::add_to_priority_queue_begin(Person* person) {
  this->priority_queue.push_front(person);
}

void Vaccine_Manager::add_to_priority_queue_end(Person* person) {
  this->priority_queue.push_back(person);
}

void Vaccine_Manager::add_to_next_dose_queue_end(Person* person) {
  this->next_dose_queue.push_back(person);
}

string Vaccine_Manager::get_vaccine_dose_priority_string() const {
  switch(this->vaccine_dose_priority) {
  case VACC_DOSE_NO_PRIORITY:
    return "No Priority";
  case VACC_DOSE_FIRST_PRIORITY:
    return "Priority, Place at Beginning of Queue";
  case VACC_DOSE_RAND_PRIORITY:
    return "Priority, Place with other Priority";
  case VACC_DOSE_LAST_PRIORITY:
    return "Priority, Place at End of Queue";
  case VACC_DOSE_SEPARATE_PRIORITY:
    return "Priority, Place at separate queue of next doses";
  default:
    FRED_WARNING("Unrecognized Vaccine Dose Priority\n");
    return "";
  }
  FRED_WARNING("Unrecognized Vaccine Dose Priority\n");
  return "";
}

void Vaccine_Manager::update(int day) {
  if(this->do_vacc) {
    this->vaccine_package->update(day);
    // Update the current vaccination capacity
    this->current_vaccine_capacity = this->vaccination_capacity_map->get_value_for_timestep(day,
											    Global::Vaccine_offset);
    cout << "Current Vaccine Stock = " << this->vaccine_package->get_vaccine(0)->get_current_stock()
	 << "\n";

    if(this->refresh_vaccine_queues_daily) {
      // update queues
      this->priority_queue.clear();
      this->queue.clear();
      fill_queues();
    }

    if(this->enable_vaccine_priority_discrete_refill == true){
      if(this->current_priority_included >= (this->policies.size() - 1)){
	if(this->vaccine_priority_timing_vector[this->current_priority_included + 1] == day || this->priority_queue.size() == 0){
	  add_next_priority_to_queues();
	}
      }
    }

    // vaccinate people in the queues:
    vaccinate(day);

    printf("After vaccinating, process immunity events\n");
    //Update events for: gaining immunity and loss of imminuty
    this->process_vaccine_next_dose_events(day);
    this->process_vaccine_immunity_start_events(day);
    this->process_vaccine_immunity_end_events(day);
  }
}

void Vaccine_Manager::reset() {
  this->priority_queue.clear();
  this->queue.clear();
  if(this->do_vacc) {
    fill_queues();
    this->vaccine_package->reset();
  }
}

void Vaccine_Manager::print() {
  this->vaccine_package->print();
}

void Vaccine_Manager::process_vaccine_immunity_start_events(int day){
  int size = this->vaccine_immunity_start_event_queue->get_size(day);
  printf("VAX_IMM_START_EVENT_QUEUE day %d size %d\n", day, size);
  if(size <= 0){
    return;
  }
  for(int i = 0; i < size; ++i) {
    Person* person =  this->vaccine_immunity_start_event_queue->get_event(day, i);

    FRED_VERBOSE(1,"vaccine_immunity_start_event day %d person %d\n",
		 day, person->get_id());
    person->update_vaccine_interventions(day);
  }
  this->vaccine_immunity_start_event_queue->clear_events(day);
}

void Vaccine_Manager::process_vaccine_next_dose_events(int day){
  int size = this->vaccine_next_dose_event_queue->get_size(day);
  printf("VAX_NEXT_DOSE_EVENT_QUEUE day %d size %d\n", day, size);
  if(size <= 0){
    return;
  }
  for(int i = 0; i < size; ++i) {
    Person* person =  this->vaccine_next_dose_event_queue->get_event(day, i);

    FRED_VERBOSE(1,"vaccine_next_dose_event day %d person %d\n",
		 day, person->get_id());
    person->update_vaccine_interventions(day);
  }
  this->vaccine_next_dose_event_queue->clear_events(day);
}

void Vaccine_Manager::process_vaccine_immunity_end_events(int day){
  int size = this->vaccine_immunity_end_event_queue->get_size(day);
  printf("VAX_IMM_END_EVENT_QUEUE day %d size %d\n", day, size);
  if(size <= 0){
    return;
  }
  for(int i = 0; i < size; ++i) {
    Person* person =  this->vaccine_immunity_end_event_queue->get_event(day, i);

    FRED_VERBOSE(1,"vaccine_immunity_end_event day %d person %d\n",
		 day, person->get_id());
    person->update_vaccine_interventions(day);
  }
  this->vaccine_immunity_end_event_queue->clear_events(day);
}

void Vaccine_Manager::vaccinate(int day) {
  if(this->do_vacc) {
    cout << "Vaccinating!\n";
  } else {
    cout << "Not vaccinating!\n";
    return;
  }

  int number_vaccinated = 0;
  int number_total_doses = 0;
  int n_p_vaccinated = 0;
  int n_r_vaccinated = 0;
  int accept_count = 0;
  int reject_count = 0;
  int reject_state_count = 0;

  // Figure out the total number of vaccines we can hand out today
  int total_vaccines_avail = this->vaccine_package->get_total_vaccines_avail_today();
  printf("Day %d Vaccine Capacity %d Total Vaccines Available %d\n", day, current_vaccine_capacity, total_vaccines_avail);
  
  if(Global::Debug > 0) {
    cout << "Vaccine Capacity on Day " << day << " = " << current_vaccine_capacity << "\n";
    cout << "Queues at beginning of vaccination:  priority (" << priority_queue.size()
	 << ")    Regular (" << this->queue.size() << ") Next dose ("<<this->next_dose_queue.size()<<")\n";
  }
  if(total_vaccines_avail == 0 || current_vaccine_capacity == 0) {
    printf("Returning. No vaccine available?(%d) on day 0 or no capacity?(%d)\n", total_vaccines_avail, current_vaccine_capacity);
    if(Global::Debug > 1) {
      cout << "No Vaccine Available on Day " << day << "\n";
    }
    Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
    Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
    Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
    Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
    Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
    return;
  }

  

  // Start vaccinating Second doses
  list<Person*>::iterator ip;
  ip = this->next_dose_queue.begin();
  while(ip != this->next_dose_queue.end()) {
    Person* current_person = *ip;
    if(current_person == NULL){
      ip = this->next_dose_queue.erase(ip);
      continue;
    }
    if(current_person->is_alive() == false){
      ip = this->next_dose_queue.erase(ip);
      continue;
    }
    if(current_person->get_health()->is_vaccinated() == 0) {
      ip = this->next_dose_queue.erase(ip);
      continue;
    }
    // Pick from vaccines of the first dose, not other vaccine!!!
    //printf("person = %d age = %.1f vacc_app = %d\n", current_person->get_id(), current_person->get_real_age(), vacc_app);
    int vacc_app = current_person->get_health()->get_vaccinated_id();
    if(vacc_app > -1){
      if(this->vaccine_package->get_vaccine(vacc_app)->get_current_stock() <= 0){
	vacc_app = -1;
      }
    }
    if(vacc_app > -1) {
      bool accept_vaccine = false;
      // STB need to refactor to work with multiple diseases
      if((this->vaccinate_symptomatics == false)
	 && (current_person->get_health()->get_symptoms_start_date(0) != -1)
	 && (day >= current_person->get_health()->get_symptoms_start_date(0))) {
	// Add an additional layer for underreporting of symptomatics
        accept_vaccine = false;
      } else {
	if(Global::Enable_Behaviors == true){
	  accept_vaccine = current_person->acceptance_of_another_vaccine_dose();	    
	}else{
	  accept_vaccine = true; // Assuming second doses are always accepted in this way
	  // printf("SECOND DOSE ACCEPTED: person = %d age = %.1f vacc_app = %d\n", current_person->get_id(), current_person->get_real_age(), vacc_app);
	}
      }
      if(accept_vaccine == true) {
        accept_count++;
        number_vaccinated++;
	//printf("PERSON accepting vaccine, accept count %d, number vaccinated %d\n", accept_count, number_vaccinated);
        this->current_vaccine_capacity--;
        n_p_vaccinated++;
        Vaccine* vacc = this->vaccine_package->get_vaccine(vacc_app);
        vacc->remove_stock(1);
        total_vaccines_avail--;
        current_person->take_vaccine(vacc, day, this);
	int curr_dose = current_person->get_current_vaccine_dose(0);
	if(curr_dose < vacc->get_number_doses() - 1 && curr_dose > -1){
	  int days_next_dose = current_person->get_days_to_next_dose(0);	  
	  if(days_next_dose > 0){
	    this->vaccine_next_dose_event_queue->add_event(days_next_dose, current_person);
	  }
	}else{
	  number_total_doses++;
	}
        ip = this->next_dose_queue.erase(ip);  // remove a vaccinated person
	// ADD VACCINE EVENTS TO QUEUE TO PROCESS LATER
	int is_vax_effective = current_person->is_vaccine_effective_any();
	if(is_vax_effective != -1){
	  int eff_day = current_person->get_vaccination_any_effective_day();
	  //printf("Day %d Vaccination is effective for person %d on day %d\n", day, current_person->get_id(), eff_day);
	  if(eff_day > -1){
	    this->vaccine_immunity_start_event_queue->add_event(eff_day, current_person);
	    
	    // Cancel immunity end events first
	    int eff_end_day = current_person->get_vaccination_immunity_loss_day();
	    if(eff_end_day > day){
	      // Add new immunity start/end events
	      // Make sure when ending immunity to check that the day is the immunity end day
	      this->vaccine_immunity_end_event_queue->add_event(eff_end_day, current_person);
	    }
	  }
	}
      } else {
        reject_count++;
        if(0) {
          ++ip;
        } else {
          // remove non-compliant person if not HBM
          ip = this->next_dose_queue.erase(ip);
        }
      }
    } else {
      if(Global::Verbose > 0) {
        cout << "Vaccine next dose not applicable for agent " << current_person->get_id() << " "
	     << current_person->get_real_age() << "\n";
      }
      ++ip;
    }

    if(total_vaccines_avail == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated next dose priority to stock out " << n_p_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  Priority (" << priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
    if(current_vaccine_capacity == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated priority to capacity " << n_p_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  Priority (" << this->priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
  }

  // If there are enough vaccines after going through the next dose queue, then:
  // Start vaccinating Priority
  //list<Person*>::iterator ip;
  ip = this->priority_queue.begin();
  // Run through the priority queue first 
  while(ip != this->priority_queue.end()) {
    Person* current_person = *ip;
    if(current_person == NULL){
      ip = this->priority_queue.erase(ip);
      continue;
    }
    if(current_person->is_alive() == false){
      ip = this->priority_queue.erase(ip);
      continue;
    }
    // If the person is in the queue, why do we need to check for age?
    int vacc_app = this->vaccine_package->pick_from_applicable_vaccines((double)(current_person->get_age()));
    //printf("person = %d age = %.1f vacc_app = %d\n", current_person->get_id(), current_person->get_real_age(), vacc_app);
    if(vacc_app > -1) {
      bool accept_vaccine = false;
      // STB need to refactor to work with multiple diseases
      if((this->vaccinate_symptomatics == false)
	 && (current_person->get_health()->get_symptoms_start_date(0) != -1)
	 && (day >= current_person->get_health()->get_symptoms_start_date(0))) {
	// Add an additional layer for underreporting of symptomatics
        accept_vaccine = false;
        reject_state_count++;
      } else {
        if(current_person->get_health()->is_vaccinated()) {
	  if(Global::Enable_Behaviors == true){
	    accept_vaccine = current_person->acceptance_of_another_vaccine_dose();	    
	  }else{
	    double r = Random::draw_random();
	    accept_vaccine = (r < this->vaccine_acceptance_prob);
	    //printf("SECOND DOSE ACCEPTED: person = %d age = %.1f vacc_app = %d\n", current_person->get_id(), current_person->get_real_age(), vacc_app);
	  }
        } else {
	  if(Global::Enable_Behaviors == true){
	    accept_vaccine = current_person->acceptance_of_vaccine();
	  }else{
	    double r = Random::draw_random();
	    accept_vaccine = (r < this->vaccine_acceptance_prob);
	  }
        }
      }
      if(accept_vaccine == true) {
        accept_count++;
        number_vaccinated++;
	//printf("PERSON accepting vaccine, accept count %d, number vaccinated %d\n", accept_count, number_vaccinated);
        this->current_vaccine_capacity--;
        n_p_vaccinated++;
        Vaccine* vacc = this->vaccine_package->get_vaccine(vacc_app);
        vacc->remove_stock(1);
        total_vaccines_avail--;
        current_person->take_vaccine(vacc, day, this);
	int curr_dose = current_person->get_current_vaccine_dose(0);
	if(curr_dose < vacc->get_number_doses() - 1 && curr_dose > -1){
	  int days_next_dose = current_person->get_days_to_next_dose(0);	  
	  if(days_next_dose > 0){
	    this->vaccine_next_dose_event_queue->add_event(days_next_dose, current_person);
	  }
	}else{
	  number_total_doses++;
	}
        ip = this->priority_queue.erase(ip);  // remove a vaccinated person
	// ADD VACCINE EVENTS TO QUEUE TO PROCESS LATER
	int is_vax_effective = current_person->is_vaccine_effective_any();
	if(is_vax_effective != -1){
	  int eff_day = current_person->get_vaccination_any_effective_day();
	  if(eff_day > -1){
	    this->vaccine_immunity_start_event_queue->add_event(eff_day, current_person);
	    
	    // Cancel immunity end events first
	    int eff_end_day = current_person->get_vaccination_immunity_loss_day();
	    if(eff_end_day > day){
	      // Add new immunity start/end events
	      // Make sure when ending immunity to check that the day is the immunity end day
	      this->vaccine_immunity_end_event_queue->add_event(eff_end_day, current_person);
	    }
	  }
	}
      } else {
        reject_count++;
	// TODO: HBM FIX THIS!
        if(0) {
          ++ip;
        } else {
          // remove non-compliant person if not HBM
          ip = this->priority_queue.erase(ip);
        }
      }
    } else {
      if(Global::Verbose > 0) {
        cout << "Vaccine not applicable for agent " << current_person->get_id() << " "
	     << current_person->get_real_age() << "\n";
      }
      ++ip;
    }

    if(total_vaccines_avail == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated priority to stock out " << n_p_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  Priority (" << priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
    if(current_vaccine_capacity == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated priority to capacity " << n_p_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  Priority (" << this->priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
  }

  if(Global::Verbose > 0) {
    cout << "Vaccinated priority to population " << n_p_vaccinated << " agents, for a total of "
	 << number_vaccinated << " on day " << day << "\n";
  }

  // Run now through the regular queue
  ip = this->queue.begin();
  while(ip != this->queue.end()) {
    Person* current_person = *ip;
    if(current_person == NULL){
      ip = this->priority_queue.erase(ip);
      continue;
    }
    if(current_person->is_alive() == false){
      ip = this->priority_queue.erase(ip);
      continue;
    }
    int vacc_app = this->vaccine_package->pick_from_applicable_vaccines(current_person->get_real_age());
    if(vacc_app > -1) {
      bool accept_vaccine = true;
      if((this->vaccinate_symptomatics == false)
	 && (current_person->get_health()->get_symptoms_start_date(0) != -1)
	 && (day >= current_person->get_health()->get_symptoms_start_date(0))) {
        accept_vaccine = false;
        reject_state_count++;
        // printf("vaccine rejected by person %d age %0.1f -- ALREADY SYMPTOMATIC\n", current_person->get_id(), current_person->get_real_age());
      } else {
        if(current_person->get_health()->is_vaccinated()) {
          // printf("vaccine rejected by person %d age %0.1f -- ALREADY VACCINATED\n", current_person->get_id(), current_person->get_real_age());
	  if(Global::Enable_Behaviors == true){
	    accept_vaccine = current_person->acceptance_of_another_vaccine_dose();
	  }else{
	    double r = Random::draw_random();
	    accept_vaccine = (r < this->vaccine_acceptance_prob);
	  }
        } else {
	  if(Global::Enable_Behaviors == true){
	    accept_vaccine = current_person->acceptance_of_vaccine();
	  }else{
	    double r = Random::draw_random();
	    accept_vaccine = (r < this->vaccine_acceptance_prob);
	  }
        }
      }
      if(accept_vaccine == true) {
        // printf("vaccine accepted by person %d age %0.1f\n", current_person->get_id(), current_person->get_real_age());
        accept_count++;
        number_vaccinated++;
        this->current_vaccine_capacity--;
        n_r_vaccinated++;	
        Vaccine* vacc = this->vaccine_package->get_vaccine(vacc_app);
        vacc->remove_stock(1);
        total_vaccines_avail--;
        current_person->take_vaccine(vacc, day, this);
	int curr_dose = current_person->get_current_vaccine_dose(0);
	if(curr_dose < vacc->get_number_doses() - 1 && curr_dose > -1){	  
	  int days_next_dose = current_person->get_days_to_next_dose(0);
	  if(days_next_dose > 0){
	    this->vaccine_next_dose_event_queue->add_event(days_next_dose, current_person);
	  }
	}else{
	  number_total_doses++;
	}
        ip = this->queue.erase(ip);  // remove a vaccinated person

	int is_vax_effective = current_person->is_vaccine_effective_any();
	if(is_vax_effective != -1){
	  int eff_day = current_person->get_vaccination_any_effective_day();
	  if(eff_day > -1){
	    this->vaccine_immunity_start_event_queue->add_event(eff_day, current_person);	    
	    // Cancel immunity end events first
	    int eff_end_day = current_person->get_vaccination_immunity_loss_day();
	    if(eff_end_day > day){
	      // Add new immunity start/end events
	      // Make sure when ending immunity to check that the day is the immunity end day
	      this->vaccine_immunity_end_event_queue->add_event(eff_end_day, current_person);
	    }
	  }
	}
      } else {
        // printf("vaccine rejected by person %d age %0.1f\n", current_person->get_id(), current_person->get_real_age());
        reject_count++;
        // skip non-compliant person under HBM
        // if(strcmp(Global::Behavior_model_type,"HBM") == 0) ip++;
        if(0)
          ip++;
        // remove non-compliant person if not HBM
        else
          ip = this->queue.erase(ip);
      }
    } else {
      ip++;
    }
    if(total_vaccines_avail == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated regular to stock_out " << n_r_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  priority (" << priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
    if(this->current_vaccine_capacity == 0) {
      if(Global::Verbose > 0) {
        cout << "Vaccinated regular to capacity " << n_r_vaccinated << " agents, for a total of "
	     << number_vaccinated << " on day " << day << "\n";
        cout << "Left in queues:  priority (" << this->priority_queue.size() << ")    Regular ("
	     << queue.size() << ")\n";
        cout << "Number of acceptances: " << accept_count << ", Number of rejections: "
	     << reject_count << "\n";
      }
      Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
      Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
      Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
      return;
    }
  }

  if(Global::Verbose > 0) {
    cout << "Vaccinated regular to population " << n_r_vaccinated << " agents, for a total of "
	 << number_vaccinated << " on day " << day << "\n";
    cout << "Left in queues:  priority (" << this->priority_queue.size() << ")    Regular ("
	 << queue.size() << ")\n";
    cout << "Number of acceptances: " << accept_count << ", Number of rejections: " << reject_count
	 << "\n";
  }
  Global::Daily_Tracker->set_index_key_pair(day,"V", number_vaccinated);
  Global::Daily_Tracker->set_index_key_pair(day,"Va", accept_count);
  Global::Daily_Tracker->set_index_key_pair(day,"Vr", reject_count);
  Global::Daily_Tracker->set_index_key_pair(day,"Vs", reject_state_count);
  Global::Daily_Tracker->set_index_key_pair(day,"Vtd", number_total_doses);
  return;
}
