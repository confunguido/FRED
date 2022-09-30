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
// File: School.h
//
#ifndef _FRED_SCHOOL_H
#define _FRED_SCHOOL_H

#include <vector>

#include "Global.h"
#include "Place.h"

#define GRADES 20

struct Time_Step_Map_Closure {
  int sim_day_start;
  int sim_day_end;
  int grade_min;
  int grade_max;
  double capacity_open;
  int income_school;
  long int census_tract;
  
  const std::string to_string() const {
    std::stringstream ss;
    ss << "School closure Time Step Map ";
    ss << " sim_day_start " << sim_day_start;
    ss << " sim_day_end " << sim_day_end;
    ss << " grade_min " << grade_min;
    ss << " grade_max " << grade_max;
    ss << " capacity_open " << capacity_open;
    ss << " income_school " << income_school;
    ss << " census_tract " << census_tract; 
    ss << std::endl;
    return ss.str();
  }
};


class Classroom;

class School : public Place {
public:
  /**
   * Default constructor
   * Note: really only used by Allocator
   */
  School();

  /**
   * Constructor with necessary parameters
   */
  School(const char* lab, char _subtype, fred::geo lon, fred::geo lat);

  ~School() {
  }

  void prepare();
  static void get_parameters();
  int get_group(int disease_id, Person* per);
  double get_transmission_prob(int disease_id, Person* i, Person* s);
  void close(int day, int day_to_close, int duration);
  void close_by_grade(int day, int day_to_close, int duration, int min_grade, int max_grade, double capacity_in);
  bool is_open(int day);
  bool should_be_open(int day, int disease_id);
  bool should_be_open_grade(int day, int grade_in);
  static bool is_global_closure_schedule_enabled(){return School::global_closure_schedule_is_enabled;}
  void apply_global_school_closure_policy(int day, int disease_id);
  void apply_individual_school_closure_policy(int day, int disease_id);
  void apply_global_schedule_school_closure_policy(int day, int disease_id);
  double get_contacts_per_day(int disease_id);
  int enroll(Person* per);
  void unenroll(int pos);
  int get_max_grade() {
    return this->max_grade;
  }

  int get_orig_students_in_grade(int grade) {
    if(grade < 0 || this->max_grade < grade) {
      return 0;
    }
    return this->orig_students_in_grade[grade];
  }

  int get_students_in_grade(int grade) {
    if(grade < 0 || this->max_grade < grade) {
      return 0;
    }
    return this->students_in_grade[grade];
  }

  int get_classrooms_in_grade(int grade) {
    if(grade < 0 || GRADES <= grade) {
      return 0;
    }
    return static_cast<int>(this->classrooms[grade].size());
  }

  void print_size_distribution();
  void print(int disease);
  int get_number_of_rooms();
  // int get_number_of_classrooms() { return (int) classrooms.size(); }
  void setup_classrooms(Allocator<Classroom> &classroom_allocator);
  Place* select_classroom_for_student(Person* per);
  int get_number_of_students() { 
    int n = 0;
    for(int grade = 1; grade < GRADES; ++grade) {
      n += this->students_in_grade[grade];
    }
    return n;
  }

  int get_orig_number_of_students() const {
    int n = 0;
    for(int grade = 1; grade < GRADES; ++grade) {
      n += this->orig_students_in_grade[grade];
    }
    return n;
  }

  static int get_max_classroom_size() {
    return School::school_classroom_size;
  }

  void set_county(int _county_index) {
    this->county_index = _county_index;
  }

  int get_county() {
    return this->county_index;
  }

  void set_income_quartile(int _income_quartile) {
    if(_income_quartile < Global::Q1 || _income_quartile > Global::Q4) {
      this->income_quartile = -1;
    } else {
      this->income_quartile = _income_quartile;
    }
  }
  
  int get_school_income(){
    return this->school_income;
  }

  int get_income_quartile() const {
    return this->income_quartile;
  }

  void set_school_income(int _school_income){
    this->school_income = _school_income;
  }

  
  void set_school_type(int _school_type){
    this->school_type = _school_type;
  }
  
  void prepare_income_quartile_pop_size() {
    if(Global::Report_Childhood_Presenteeism) {
      int size = get_size();
      // update population stats based on income quartile of this school
      if(this->income_quartile == Global::Q1) {
        School::pop_income_Q1 += size;
      } else if(this->income_quartile == Global::Q2) {
        School::pop_income_Q2 += size;
      } else if(this->income_quartile == Global::Q3) {
        School::pop_income_Q3 += size;
      } else if(this->income_quartile == Global::Q4) {
        School::pop_income_Q4 += size;
      }
      School::total_school_pop += size;
    }
  }

  static int get_total_school_pop() {
    return School::total_school_pop;
  }

  static int get_school_pop_income_quartile_1() {
    return School::pop_income_Q1;
  }

  static int get_school_pop_income_quartile_2() {
    return School::pop_income_Q2;
  }

  static int get_school_pop_income_quartile_3() {
    return School::pop_income_Q3;
  }

  static int get_school_pop_income_quartile_4() {
    return School::pop_income_Q4;
  }

  static char* get_school_closure_policy() {
    return School::school_closure_policy;
  }

  //for access from Classroom:
  static double get_school_contacts_per_day(int disease_id) {
    return School::contacts_per_day;
  }

private:
  static double contacts_per_day;
  static double** prob_transmission_per_contact;
  static char school_closure_policy[];
  static int school_closure_day;
  static int min_school_closure_day;
  static double school_closure_threshold;
  static double individual_school_closure_threshold;
  static double individual_school_closure_current_case_rate_threshold;
  static int school_closure_cases;
  static int school_closure_current_cases;
  static int school_closure_duration;
  static int school_closure_delay;
  static bool individual_school_closure_by_cases;
  static bool individual_school_closure_by_current_cases;
  static bool individual_school_closure_by_wastewater;
  static bool school_include_rna_measurement_variability;
  static int individual_school_wastewater_threshold;
  static double school_wastewater_measurement_negbin_size;
  static double school_underreporting_rate;
  static int school_wastewater_measurement_frequency;
  static int school_case_reporting_frequency;
  static int school_summer_schedule;
  static char school_summer_start[];
  static char school_summer_end[];
  static int summer_start_month;
  static int summer_start_day;
  static int summer_end_month;
  static int summer_end_day;
  static int school_classroom_size;
  static int school_classroom_size_array[GRADES];
  static bool global_closure_is_active;
  static bool global_closure_schedule_is_enabled;
  static int global_close_date;
  static int global_open_date;

  static int total_school_pop;
  static int pop_income_Q1;
  static int pop_income_Q2;
  static int pop_income_Q3;
  static int pop_income_Q4;
  static std::vector<Time_Step_Map_Closure * > school_closure_schedule;
  int students_in_grade[GRADES];

  int orig_students_in_grade[GRADES];
  int next_classroom[GRADES];
  vector<Classroom*> classrooms[GRADES];
  bool closure_dates_have_been_set;
  bool closure_grade_dates_have_been_set;
  bool school_closure_dates_checked_today;
  int day_closure_dates_grades_set;
  int max_grade;
  int county_index;
  int income_quartile;
  int wastewater_rna;
  int wastewater_rna_actual;
  int school_current_symptomatic_infections;
  int last_wastewater_measurement_day;
  int last_case_reporting_day;
  int last_current_case_reporting_day;
  int last_current_symptomatic_infections_update_day;
  int close_grade_date[GRADES];
  int open_grade_date[GRADES];
  double open_capacity_grade[GRADES];
  int school_income;
  int school_type;
};

#endif // _FRED_SCHOOL_H
