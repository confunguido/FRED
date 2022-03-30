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
// File: School.cc
//
#include <cmath>
#include "School.h"
#include "Classroom.h"
#include "Date.h"
#include "Disease.h"
#include "Disease_List.h"
#include "Global.h"
#include "Params.h"
#include "Person.h"
#include "Place_List.h"
#include "Population.h"
#include "Random.h"
#include "Utils.h"
#include "Tracker.h"

//Private static variables that will be set by parameter lookups
double School::contacts_per_day;
double** School::prob_transmission_per_contact;
int School::school_classroom_size = 0;
char School::school_closure_policy[80] = "undefined";
int School::school_closure_day = 0;
int School::min_school_closure_day = 0;
double School::school_closure_threshold = 0.0;
double School::individual_school_closure_threshold = 0.0;
int School::school_closure_cases = -1;
int School::school_closure_duration = 0;
int School::school_closure_delay = 0;
bool School::individual_school_closure_by_cases = false;
bool School::individual_school_closure_by_wastewater = false;
bool School::school_include_rna_measurement_variability = false;
int School::individual_school_wastewater_threshold = 0;
double School::school_wastewater_measurement_negbin_size = INFINITY;
int School::school_summer_schedule = 0;
char School::school_summer_start[8];
char School::school_summer_end[8];
int School::summer_start_month = 0;
int School::summer_start_day = 0;
int School::summer_end_month = 0;
int School::summer_end_day = 0;
int School::total_school_pop = 0;
int School::pop_income_Q1 = 0;
int School::pop_income_Q2 = 0;
int School::pop_income_Q3 = 0;
int School::pop_income_Q4 = 0;
bool School::global_closure_is_active = false;
int School::global_close_date = 0;
int School::global_open_date = 0;
bool School::global_closure_schedule_is_enabled = false;

int School::school_classroom_size_array[GRADES];

std::vector<Time_Step_Map_Closure * > School::school_closure_schedule;

School::School() : Place() {
  this->set_type(Place::TYPE_SCHOOL);
  this->set_subtype(Place::SUBTYPE_NONE);
  this->intimacy = 0.025;
  for(int i = 0; i < GRADES; ++i) {
    this->students_in_grade[i] = 0;
    this->orig_students_in_grade[i] = 0;
    this->next_classroom[i] = 0;
    this->classrooms[i].clear();
    this->close_grade_date[i] = INT_MAX;
    this->open_grade_date[i] = 0;
    this->open_capacity_grade[i] = 1.0;
  }
  this->closure_dates_have_been_set = false;
  this->closure_grade_dates_have_been_set = false;
  this->school_closure_dates_checked_today = false;
  this->day_closure_dates_grades_set = -1;
  this->staff_size = 0;
  this->max_grade = -1;
  this->county_index = -1;
  this->income_quartile = -1;
  this->school_income = 0;
  this->school_type = 0;
}


School::School(const char* lab, char _subtype, fred::geo lon, fred::geo lat) : Place(lab, lon, lat) {
  this->set_type(Place::TYPE_SCHOOL);
  this->set_subtype(_subtype);
  this->intimacy = 0.025;
  for(int i = 0; i < GRADES; ++i) {
    this->students_in_grade[i] = 0;
    this->orig_students_in_grade[i] = 0;
    this->next_classroom[i] = 0;
    this->classrooms[i].clear();
  }
  this->closure_dates_have_been_set = false;
  this->closure_grade_dates_have_been_set = false;
  this->school_closure_dates_checked_today = false;
  this->day_closure_dates_grades_set = -1;
  this->staff_size = 0;
  this->max_grade = -1;
  this->county_index = -1;
  this->income_quartile = -1;
  this->school_income = 0;
  this->school_type = 0;
}

void School::prepare() {
  assert(Global::Pop.is_load_completed());
  // call base class function to perform preparations common to all Places 
  Place::prepare();

  // record original size in each grade
  for(int g = 0; g < GRADES; ++g) {
    this->orig_students_in_grade[g] = this->students_in_grade[g];
  }
}


void School::get_parameters() {

  Params::get_param_from_string("school_contacts", &School::contacts_per_day);
  int n = Params::get_param_matrix((char *)"school_trans_per_contact", &School::prob_transmission_per_contact);
  if(Global::Verbose > 1) {
    printf("\nSchool_contact_prob:\n");
    for(int i  = 0; i < n; ++i)  {
      for(int j  = 0; j < n; ++j) {
	      printf("%f ", School::prob_transmission_per_contact[i][j]);
      }
      printf("\n");
    }
  }

  // normalize contact parameters
  // find max contact prob
  double max_prob = 0.0;
  for(int i  = 0; i < n; ++i)  {
    for(int j  = 0; j < n; ++j) {
      if (School::prob_transmission_per_contact[i][j] > max_prob) {
	      max_prob = School::prob_transmission_per_contact[i][j];
      }
    }
  }

  // convert max contact prob to 1.0
  if(max_prob > 0) {
    for(int i  = 0; i < n; ++i)  {
      for(int j  = 0; j < n; ++j) {
	      School::prob_transmission_per_contact[i][j] /= max_prob;
      }
    }
    // compensate contact rate
    School::contacts_per_day *= max_prob;
  }

  if(Global::Verbose > 0) {
    printf("\nSchool_contact_prob after normalization:\n");
    for(int i  = 0; i < n; ++i)  {
      for(int j  = 0; j < n; ++j) {
	      printf("%f ", School::prob_transmission_per_contact[i][j]);
      }
      printf("\n");
    }
    printf("\ncontact rate: %f\n", School::contacts_per_day);
  }
  // end normalization
  
  Params::get_param_from_string("school_classroom_size", &School::school_classroom_size);

  // Specific classroom sizes for each grade
  if(Global::Enable_School_Classroom_Size_Array == true){
    char class_size_str[MAX_PARAM_SIZE];
    std::vector<double> classroom_size_arr;
    Params::get_param((char*)"school_classroom_size_array",class_size_str);
    Params::get_param_vector_from_string(class_size_str,classroom_size_arr);
    if(classroom_size_arr.size() != GRADES){
      Utils::fred_abort("Help!  classroom sizes are not properly setup %s \n", class_size_str);
      abort();
    }
    for(int cc = 0; cc < GRADES; cc++){
      School::school_classroom_size_array[cc] = (int) classroom_size_arr[cc];
    }
    for(int cc = 0; cc < GRADES; cc++){
      printf("CLASSROOM FOR GRADE %d has size %d\n", cc, school_classroom_size_array[cc]);
    }
  }
    
  // summer school parameters
  Params::get_param_from_string("school_summer_schedule", &School::school_summer_schedule);
  Params::get_param_from_string("school_summer_start", School::school_summer_start);
  Params::get_param_from_string("school_summer_end", School::school_summer_end);
  sscanf(School::school_summer_start, "%d-%d", &School::summer_start_month, &School::summer_start_day);
  sscanf(School::school_summer_end, "%d-%d", &School::summer_end_month, &School::summer_end_day);

  // school closure parameters
  Params::get_param_from_string("school_closure_policy", School::school_closure_policy);
  Params::get_param_from_string("school_closure_duration", &School::school_closure_duration);
  Params::get_param_from_string("school_closure_delay", &School::school_closure_delay);
  Params::get_param_from_string("school_closure_day", &School::school_closure_day);
  Params::get_param_from_string("min_school_closure_day", &School::min_school_closure_day);
  Params::get_param_from_string("school_closure_ar_threshold", &School::school_closure_threshold);
  Params::get_param_from_string("individual_school_closure_ar_threshold",
				&School::individual_school_closure_threshold);
  Params::get_param_from_string("school_closure_cases", &School::school_closure_cases);
  int temp_integer = 0;
  Params::get_param_from_string("individual_school_closure_by_cases",
				&temp_integer);
  School::individual_school_closure_by_cases = (temp_integer == 0 ? false : true);
  Params::get_param_from_string("individual_school_closure_by_wastewater",
				&temp_integer);
  School::individual_school_closure_by_wastewater = (temp_integer == 0 ? false : true);
  Params::get_param_from_string("school_include_rna_measurement_variability",
				&temp_integer);
  School::school_include_rna_measurement_variability = (temp_integer == 0 ? false : true);
  Params::get_param_from_string("individual_school_wastewater_threshold",
				&School::individual_school_wastewater_threshold);
  Params::get_param_from_string("school_wastewater_measurement_negbin_size",
				&School::school_wastewater_measurement_negbin_size);

  // if global_schedule_file is not none, then read it
  if(strcmp(School::school_closure_policy, "global_schedule") == 0) {
    School::global_closure_schedule_is_enabled = true;
    char map_file_name[FRED_STRING_SIZE];
    Params::get_param_from_string("school_global_schedule_file", map_file_name);
    printf("READING SCHOOL SCHEDULE %s\n", map_file_name);
    // If this parameter is "none", then there is no map
    if(strncmp(map_file_name, "none", 4) != 0){
      Utils::get_fred_file_name(map_file_name);
      printf("READING: School schedule file: %s\n", map_file_name);
      ifstream* ts_input = new ifstream(map_file_name);
      if(!ts_input->is_open()) {
	Utils::fred_abort("Help!  Can't read %s School Timestep Map\n", map_file_name);
	abort();
      }
      string line;
      while(getline(*ts_input,line)){		
	if(line[0] == ' ' || line[0] == '\n' || line[0] == '#') { // empty line or comment
	  continue;
	}
	char cstr[FRED_STRING_SIZE];
	std::strcpy(cstr, line.c_str());
	Time_Step_Map_Closure * tmap = new Time_Step_Map_Closure;
	int n = sscanf(cstr,
		       "%d %d %d %d %lg %d %ld",
		       &tmap->sim_day_start, &tmap->sim_day_end, &tmap->grade_min, &tmap->grade_max, &tmap->capacity_open, &tmap->income_school, &tmap->census_tract);
	//printf("SCHOOL SCHEDULE LINES: %d\n",n);
	if(n < 5) {
	  Utils::fred_abort("Need to specify at least SimulationDayStart, SimulationDayEnd, Minimum grade, and Maximum Grade (0-20), capacity open");
	}
	if(n < 7){
	  tmap->census_tract = 0;
	  if(n < 6){
	    tmap->income_school = 0;
	  }
	}

	if(tmap->capacity_open < 0.0){tmap->capacity_open = 0.0;}
	if(tmap->capacity_open > 1.0){tmap->capacity_open = 1.0;}
	if(tmap->sim_day_end > tmap->sim_day_start && tmap->grade_max > tmap->grade_min){
	  if(tmap->grade_min > 0 && tmap->grade_min <= GRADES){
	    if(tmap->grade_max > 0 && tmap->grade_max <= GRADES){
	      if(tmap->grade_max > tmap->grade_min){
		School::school_closure_schedule.push_back(tmap);
	      }
	    }
	  }
	}
      }
      ts_input->close();
    }
  }
  /*
  for(int i = 0; i < School::school_closure_schedule.size(); ++i){
    string ss = School::school_closure_schedule[i]->to_string();
    printf("%s\n", ss.c_str());
    }*/
    
  
  
  // aliases for parameters
  int Weeks;
  Params::get_param_from_string("Weeks", &Weeks);
  if(Weeks > -1) {
    School::school_closure_duration = 7 * Weeks;
  }

  int cases;
  Params::get_param_from_string("Cases", &cases);
  if(cases > -1) {
    School::school_closure_cases = cases;
  }
}

int School::get_group(int disease_id, Person* per) {
  int age = per->get_age();
  if(age < 12) {
    return 0;
  } else if(age < 16) {
    return 1;
  } else if(per->is_student()) {
    return 2;
  } else {
    return 3;
  }
}

double School::get_transmission_prob(int disease_id, Person* i, Person* s) {
  // i = infected agent
  // s = susceptible agent
  int row = get_group(disease_id, i);
  int col = get_group(disease_id, s);
  double tr_pr = School::prob_transmission_per_contact[row][col];
  return tr_pr;
}

void School::close(int day, int day_to_close, int duration) {
  this->close_date = day_to_close;
  this->open_date = close_date + duration;
  this->closure_dates_have_been_set = true;
  
  // log this school closure decision
  if(Global::Verbose > 0) {
    printf("SCHOOL %s CLOSURE decision day %d close_date %d duration %d open_date %d\n",
	   this->get_label(), day, this->close_date, duration, this->open_date);
  }
}

void School::close_by_grade(int day, int day_to_close, int duration, int min_grade, int max_grade, double capacity_in) {  
  for(int i = min_grade - 1; i < max_grade; ++i){
    this->close_grade_date[i] = day_to_close;
    this->open_grade_date[i] = day_to_close + duration;
    this->open_capacity_grade[i] = capacity_in;
    this->closure_grade_dates_have_been_set = true;
  }
  // log this school closure decision
  if(Global::Verbose > 0) {
    printf("SCHOOL %s CLOSURE GRADE decision day %d close_date %d duration %d open_date %d grades: %d-%d capacity %.2f\n",
	   this->get_label(), day,day_to_close, duration, close_grade_date + duration, min_grade, max_grade, capacity_in);
  }
}


bool School::is_open(int day) {
  // Ignore closure if school is open by grades
  int sch_census_tract_index = this->get_census_tract_index();
  long int sch_census_tract = -1;
  if(sch_census_tract_index > 0){
    sch_census_tract = Global::Places.get_census_tract_with_index(sch_census_tract_index);
  }
  /*if(sch_census_tract == 11001008526){
    printf("SCHOOL %s CENSUS TRACT %lu verifying if it is open on day %d\n", this->get_label(),sch_census_tract, day);
    }*/
  if(this->closure_grade_dates_have_been_set == true){
    bool open_grade = false;
    for(int grade = 0; grade < GRADES; grade++){
      if(day >= this->close_grade_date[grade] && this->open_grade_date[grade] >= day){
	open_grade = true;
	/*
	if(sch_census_tract == 11001008526){
	  printf("SCHOOL CENSUS TRACT %lu open for grade %d\n", sch_census_tract, grade);
	  }*/
	break;
      }
    }
    if(open_grade == true){
      /*
      if(sch_census_tract == 11001008526){
      printf("SCHOOL CENSUS TRACT %lu open at least for one grade on day %d\n", sch_census_tract, day);
	}*/
      //printf("SCHOOL CENSUS TRACT %lu open at least for one grade on day %d\n", sch_census_tract, day);
      return true;
    }
  }
  
  bool open = (day < this->close_date || this->open_date <= day);
  if(this->open_date <= day) {
    // if(open == true) {
    this->closure_dates_have_been_set = false;
  }
  if(!open) {
    /*
    if(sch_census_tract == 11001008526){
      printf("SCHOOL CENSUS TRACT %lu closed on day %d\n", sch_census_tract, day);
      }*/
    FRED_VERBOSE(1, "Place %s is closed on day %d\n", this->get_label(), day);
  }
  /*
  else{
    if(sch_census_tract == 11001008526){
    printf("SCHOOL CENSUS TRACT %lu open for ALL GRADES on day %d, close date (%d) open date(%d)\n", sch_census_tract, day, this->close_date, this->open_date);
    }
  }
  */
  if(open){
    //printf("SCHOOL CENSUS TRACT %lu open for ALL GRADES on day %d, close date (%d) open date(%d)\n", sch_census_tract, day, this->close_date, this->open_date);
  }
  return open;
}

bool School::should_be_open_grade(int day, int grade) {
  if(School::global_closure_schedule_is_enabled == false){
    return true;
  }
  //printf("Should school %s be open day %d grade %d\n", this->get_label(), day, grade);
  
  if(this->open_capacity_grade[grade] == 1.0){
    return true;
  }
  bool open = (day < this->close_grade_date[grade] || this->open_grade_date[grade] <= day);
  if(!open){
    if(Random::draw_random() < this->open_capacity_grade[grade]) {
      open = true;
    }
  }
  if(!open) {
    FRED_VERBOSE(2, "School Place %s is closed on day %d for grade %d\n", this->get_label(), day, grade);
  }
  return open;
}

bool School::should_be_open(int day, int disease_id) {
  // no students
  if(get_size() == 0) {
    return false;
  }

  // summer break
  if(School::school_summer_schedule > 0) {
    int month = Date::get_month();
    int day_of_month = Date::get_day_of_month();
    if((month == School::summer_start_month && School::summer_start_day <= day_of_month) ||
       (month == School::summer_end_month && day_of_month <= School::summer_end_day) ||
       (School::summer_start_month < month && month < School::summer_end_month)) {
      if(Global::Verbose > 1) {
	      fprintf(Global::Statusfp, "School %s closed for summer\n", this->get_label());
	      fflush(Global::Statusfp);
      }
      return false;
    }
  }
  
  // global school closure policy in effect
  if(strcmp(School::school_closure_policy, "global_schedule") == 0) {
    if(day > this->day_closure_dates_grades_set){
      apply_global_schedule_school_closure_policy(day, disease_id);
    }    
    return is_open(day);
  }
  
  // stick to previously made decision to close
  if(this->closure_dates_have_been_set) {
    return is_open(day);
  }

  // global school closure policy in effect
  if(strcmp(School::school_closure_policy, "global") == 0) {
    apply_global_school_closure_policy(day, disease_id);
    return is_open(day);
  }

  // individual school closure policy in effect
  if(strcmp(School::school_closure_policy, "individual") == 0) {
    apply_individual_school_closure_policy(day, disease_id);
    return is_open(day);
  }  

  // if school_closure_policy is not recognized, then open
  return true;
}

void School::apply_global_schedule_school_closure_policy(int day, int disease_id) {
  // Should I move this to place list?
  // Every day, check is day is within school closure of time step
  this->closure_grade_dates_have_been_set = false;
  for(int i = 0; i < School::school_closure_schedule.size(); ++i) {
    Time_Step_Map_Closure* tmap = School::school_closure_schedule[i];
    if(tmap->sim_day_start <= day && day <= tmap->sim_day_end) {
      // if min grade = 0 and max grade = 20, then close the school, else close some grades only
      if(tmap->grade_min <= 1 && tmap->grade_max == GRADES && tmap->capacity_open == 0.0 && tmap->income_school < 1 && tmap->census_tract < 1){
	int sch_census_tract_index = this->get_census_tract_index();
	long int sch_census_tract = -1;
	if(sch_census_tract_index > 0){
	  sch_census_tract = Global::Places.get_census_tract_with_index(sch_census_tract_index);
	}
	if(sch_census_tract == 11001008526){
	  printf("Day: %d School in census_tract %ld closed with tmap %s\n", day, sch_census_tract, tmap->to_string().c_str());
	}
	close(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start);
      }else if(tmap->grade_min <= 1 && tmap->grade_max == GRADES && tmap->capacity_open == 0.0 && tmap->income_school > 0){
	if(tmap->income_school == this->school_income){
	  close(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start);
	}
      }else if(tmap->grade_min <= 1 && tmap->grade_max == GRADES && tmap->capacity_open == 0.0 && tmap->census_tract > 0){
	int sch_census_tract_index = this->get_census_tract_index();
	if(sch_census_tract_index > 0){
	  long int sch_census_tract = Global::Places.get_census_tract_with_index(sch_census_tract_index);
	  printf("Applying GLOBAL SCHEDULE with census tract: %lu, checking with school %lu\n", tmap->census_tract, sch_census_tract);
	  if(tmap->census_tract == sch_census_tract){
	    close(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start);
	  }
	}
      }else{
	if(tmap->income_school > 0){
	  if(tmap->income_school == this->school_income){
	    close_by_grade(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start, tmap->grade_min, tmap->grade_max, tmap->capacity_open);
	  }
	}else if(tmap->census_tract > 0){
	  int sch_census_tract_index = this->get_census_tract_index();
	  if(sch_census_tract_index > 0){
	    long int sch_census_tract = Global::Places.get_census_tract_with_index(sch_census_tract_index);
	    if(tmap->census_tract == sch_census_tract){
	      if(sch_census_tract == 11001008526){
		printf("Day: %d School in census_tract %ld closed with tmap %s\n", day, sch_census_tract, tmap->to_string().c_str());
	      }
	      close_by_grade(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start, tmap->grade_min, tmap->grade_max, tmap->capacity_open);
	    }
	  }
	}else{
	  close_by_grade(day,tmap->sim_day_start, tmap->sim_day_end - tmap->sim_day_start, tmap->grade_min, tmap->grade_max, tmap->capacity_open);
	}
      }
      if(Global::Verbose > 0) {
	printf("GLOBAL SCHEDULE SCHOOL CLOSURE %s\n", tmap->to_string().c_str());
      }
    }
  }
  this->day_closure_dates_grades_set = day;
  this->school_closure_dates_checked_today = false;
}

void School::apply_global_school_closure_policy(int day, int disease_id) {

  // Only test triggers for school closure if not global closure is not already activated
  if(School::global_closure_is_active == false) {
    // Setting school_closure_day > -1 overrides other global triggers.
    // Close schools if the closure date has arrived (after a delay)
    if(School::school_closure_day > -1) {
      if(School::school_closure_day <= day) {
        // the following only happens once
        School::global_close_date = day + School::school_closure_delay;
        School::global_open_date = day + School::school_closure_delay
	        + School::school_closure_duration;
	      School::global_closure_is_active = true;
      }
    } else {
      // Close schools if the global symptomatic attack rate has reached the threshold (after a delay)
      Disease* disease = Global::Diseases.get_disease(disease_id);
      if(School::school_closure_threshold <= disease->get_symptomatic_attack_rate()) {
        // the following only happens once
        School::global_close_date = day + School::school_closure_delay;
        School::global_open_date = day + School::school_closure_delay
	        + School::school_closure_duration;
	      School::global_closure_is_active = true;
      }
    }
  }
  if(School::global_closure_is_active) {
    // set close and open dates for this school (only once)
    close(day,School::global_close_date, School::school_closure_duration);

    // log this school closure decision
    if(Global::Verbose > 0) {
      Disease* disease = Global::Diseases.get_disease(disease_id);
      printf("GLOBAL SCHOOL CLOSURE pop_ar %5.2f local_cases = %d / %d (%5.2f)\n",
	     disease->get_symptomatic_attack_rate(), get_total_cases(disease_id),
	     get_size(), get_symptomatic_attack_rate(disease_id));
    }
  }
}

void School::apply_individual_school_closure_policy(int day, int disease_id) {

  // don't apply any policy prior to School::min_school_closure_day
  if(day <= School::min_school_closure_day) {
    return;
  }

  // don't apply any policy before the epidemic reaches a noticeable threshold
  Disease* disease = Global::Diseases.get_disease(disease_id);
  if(disease->get_symptomatic_attack_rate() < School::school_closure_threshold) {
    return;
  }

  bool close_this_school = false;
  bool false_negative = false;
  
  // If individual_school_closure_by_cases then close if the relevant threshold is met
  // if school_closure_cases > -1 then close if this number of cases occurs
  if(School::individual_school_closure_by_cases){
    if(School::school_closure_cases != -1) {
      close_this_school = (School::school_closure_cases <= get_total_cases(disease_id));
    } else {
      // close if attack rate threshold is reached
      close_this_school = (School::individual_school_closure_threshold <= get_symptomatic_attack_rate(disease_id));
    }
  }
  
  // double wastewater_rna;
  // if individual_school_closure_by_wastewater then close if this threshold is met
  if(School::individual_school_closure_by_wastewater) {
    wastewater_rna_actual = get_wastewater_rna(disease_id,day);
    wastewater_rna = (school_include_rna_measurement_variability ?
		      get_wastewater_rna(disease_id,day,school_wastewater_measurement_negbin_size) :
		      wastewater_rna_actual);
    close_this_school = (School::individual_school_wastewater_threshold <= wastewater_rna);
    false_negative = (!close_this_school && get_current_infections(day,disease_id) > 0);

    if (get_number_of_infectious_people(disease_id) > get_current_infections(day,disease_id)) {
      printf("ERROR!! School has more infectious people than infections...");
    }
  }

  if(close_this_school) {
    // set close and open dates for this school (only once)
    close(day,day + School::school_closure_delay, School::school_closure_duration);

    // log this school closure decision
    if(Global::Verbose > 0) {
      Disease* disease = Global::Diseases.get_disease(disease_id);
      if(School::individual_school_closure_by_wastewater) {
	printf("LOCAL SCHOOL CLOSURE pop_ar %.3f local_cases = %d / %d (%.3f) WW conc: measurement = %d GC/l, actual = %d GC/l infections = %d infectious people = %d\n",
	       disease->get_symptomatic_attack_rate(), get_total_cases(disease_id),
	       get_size(), get_symptomatic_attack_rate(disease_id),
	       wastewater_rna,wastewater_rna_actual,get_current_infections(day,disease_id),
	       get_number_of_infectious_people(disease_id));
      } else {
	printf("LOCAL SCHOOL CLOSURE pop_ar %.3f local_cases = %d / %d (%.3f)\n",
	       disease->get_symptomatic_attack_rate(), get_total_cases(disease_id),
	       get_size(), get_symptomatic_attack_rate(disease_id));
      }
    }
  } else if (false_negative && Global::Verbose > 0) {
    printf("FALSE NEGATIVE LOCAL SCHOOL CLOSURE WW conc: measurement = %d GC/l, actual = %d GC/l infections = %d infectious people = %d\n",
	   wastewater_rna,wastewater_rna_actual,get_current_infections(day,disease_id),
	   get_number_of_infectious_people(disease_id));
  }
}

double School::get_contacts_per_day(int disease_id) {
  return School::contacts_per_day;
}

int School::enroll(Person* person) {

  // call base class method:
  int return_value = Place::enroll(person);

  FRED_VERBOSE(1,"Enroll person %d age %d in school %d %s new size %d\n",
	       person->get_id(), person->get_age(), this->get_id(), this->get_label(), this->get_size());
  if(person->is_teacher()) {
    this->staff_size++;
  } else {
    int age = person->get_age();
    if (age >= GRADES){
      age == GRADES -1;
    }
    int grade = ((age < GRADES) ? age : GRADES - 1);
    assert(grade > 0);
    this->students_in_grade[grade]++;
    if(grade > this->max_grade) {
      this->max_grade = grade;
    }
    person->set_grade(grade);
  }

  return return_value;
}

void School::unenroll(int pos) {
  int size = this->enrollees.size();
  assert(0 <= pos && pos < size);
  Person* removed = this->enrollees[pos];
  int grade = removed->get_grade();
  FRED_VERBOSE(1,"UNENROLL person %d age %d grade %d, is_teacher %d from school %d %s Size = %d\n",
	       removed->get_id(), removed->get_age(), grade, removed->is_teacher()?1:0, this->get_id(), this->get_label(), this->get_size());

  // call the base class method
  Place::unenroll(pos);

  if(removed->is_teacher() || grade == 0) {
    this->staff_size--;
  } else {
    assert(0 < grade && grade <= this->max_grade);
    this->students_in_grade[grade]--;
  }
  removed->set_grade(0);
  FRED_VERBOSE(1,"UNENROLLED from school %d %s size = %d\n", this->get_id(), this->get_label(), this->get_size());
}

void School::print(int disease_id) {
  fprintf(Global::Statusfp, "Place %d label %s type %c ", this->get_id(), this->get_label(), this->get_type());
  for(int g = 0; g < GRADES; ++g) {
    fprintf(Global::Statusfp, "%d students in grade %d | ", this->students_in_grade[g], g);
  }
  fprintf(Global::Statusfp, "\n");
  fflush(Global::Statusfp);
}

int School::get_number_of_rooms() {
  int total_rooms = 0;

  for(int a = 0; a < GRADES; ++a) {
    int class_size_tmp = School::school_classroom_size;
    if(Global::Enable_School_Classroom_Size_Array){
      class_size_tmp = School::school_classroom_size_array[a];
    }
    if(class_size_tmp == 0){
      continue;
    }
    int n = this->students_in_grade[a];
    if(n == 0) {
      continue;
    }
    int rooms = n / class_size_tmp;
    if(n % class_size_tmp) {
      rooms++;
    }
    total_rooms += rooms;
  }
  return total_rooms;
}

void School::setup_classrooms(Allocator<Classroom> &classroom_allocator) {
  for(int a = 0; a < GRADES; ++a) {
    int class_size_tmp = School::school_classroom_size;
    if(Global::Enable_School_Classroom_Size_Array){
      class_size_tmp = School::school_classroom_size_array[a];
    }
    if(class_size_tmp == 0){
      continue;
    }
    int n = this->students_in_grade[a];
    if(n == 0) {
      continue;
    }
    int rooms = n / class_size_tmp;
    if(n % class_size_tmp) {
      rooms++;
    }

    FRED_STATUS(1, "School %d %s grade %d number %d rooms %d\n", this->get_id(), this->get_label(), a, n, rooms);

    for(int c = 0; c < rooms; ++c) {
      char new_label[128];
      sprintf(new_label, "%s-%02d-%02d", this->get_label(), a, c + 1);

      Classroom* clsrm = new (classroom_allocator.get_free()) Classroom(new_label,
									Place::SUBTYPE_NONE,
									this->get_longitude(),
									this->get_latitude());
      clsrm->set_school(this);

      this->classrooms[a].push_back(clsrm);
    }
  }
}

Place* School::select_classroom_for_student(Person* per) {
  if(School::school_classroom_size == 0 && Global::Enable_School_Classroom_Size_Array == false) {
    return NULL;
  }
  int grade = per->get_age();
  if(GRADES <= grade) {
    grade = GRADES - 1;
  }
  if(Global::Verbose > 1) {
    fprintf(Global::Statusfp, "assign classroom for student %d in school %d %s grade %d\n",
	    per->get_id(), this->get_id(), this->get_label(), grade);
    fflush(Global::Statusfp);
  }

  if(this->classrooms[grade].size() == 0) {
    return NULL;
  }

  // pick next classroom for this grade, round-robin
  int room = this->next_classroom[grade];
  if(room < (int) this->classrooms[grade].size() - 1) {
    this->next_classroom[grade]++;
  } else {
    this->next_classroom[grade] = 0;
  }

  // pick next classroom for this grade at random
  // int room = Random::draw_random_int(0,(classrooms[grade].size()-1));

  if(Global::Verbose > 1) {
    fprintf(Global::Statusfp, "room = %d %s %d\n", room, this->classrooms[grade][room]->get_label(),
	    this->classrooms[grade][room]->get_id());
    fflush(Global::Statusfp);
  }
  return this->classrooms[grade][room];
}

void School::print_size_distribution() {
  for(int g = 1; g < GRADES; ++g) {
    if(this->orig_students_in_grade[g] > 0) {
      printf("SCHOOL %s grade %d orig %d current %d\n", 
	     this->get_label(), g, this->orig_students_in_grade[g], this->students_in_grade[g]);
    }
  }
  fflush(stdout);
}
