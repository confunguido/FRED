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
// File: Health.cc
//
#include <new>
#include <stdexcept>

#include "Age_Map.h"
#include "Antiviral.h"
#include "AV_Health.h"
#include "AV_Manager.h"
#include "Date.h"
#include "Disease.h"
#include "Disease_List.h"
#include "Evolution.h"
#include "Infection.h"
#include "Health.h"
#include "Household.h"
#include "Manager.h"
#include "Mixing_Group.h"
#include "Past_Infection.h"
#include "Person.h"
#include "Place.h"
#include "Place_List.h"
#include "Population.h"
#include "Random.h"
#include "Utils.h"
#include "Vaccine.h"
#include "Vaccine_Dose.h"
#include "Vaccine_Health.h"
#include "Vaccine_Manager.h"

// static variables
int Health::nantivirals = -1;
Age_Map* Health::asthma_prob = NULL;
Age_Map* Health::COPD_prob = NULL;
Age_Map* Health::chronic_renal_disease_prob = NULL;
Age_Map* Health::diabetes_prob = NULL;
Age_Map* Health::heart_disease_prob = NULL;
Age_Map* Health::hypertension_prob = NULL;
Age_Map* Health::hypercholestrolemia_prob = NULL;

Age_Map* Health::asthma_hospitalization_prob_mult = NULL;
Age_Map* Health::COPD_hospitalization_prob_mult = NULL;
Age_Map* Health::chronic_renal_disease_hospitalization_prob_mult = NULL;
Age_Map* Health::diabetes_hospitalization_prob_mult = NULL;
Age_Map* Health::heart_disease_hospitalization_prob_mult = NULL;
Age_Map* Health::hypertension_hospitalization_prob_mult = NULL;
Age_Map* Health::hypercholestrolemia_hospitalization_prob_mult = NULL;

Age_Map* Health::asthma_case_fatality_prob_mult = NULL;
Age_Map* Health::COPD_case_fatality_prob_mult = NULL;
Age_Map* Health::chronic_renal_disease_case_fatality_prob_mult = NULL;
Age_Map* Health::diabetes_case_fatality_prob_mult = NULL;
Age_Map* Health::heart_disease_case_fatality_prob_mult = NULL;
Age_Map* Health::hypertension_case_fatality_prob_mult = NULL;
Age_Map* Health::hypercholestrolemia_case_fatality_prob_mult = NULL;

Age_Map* Health::pregnancy_hospitalization_prob_mult = NULL;
Age_Map* Health::pregnancy_case_fatality_prob_mult = NULL;

bool Health::is_initialized = false;

// health protective behavior parameters
int Health::Min_age_wear_face_masks = 0;
int Health::Days_to_wear_face_masks = 0;
int Health::Day_start_wearing_face_masks = 0;
std::unordered_map<string,double> Health::Face_mask_compliance;
double Health::Hand_washing_compliance = 0.0;

double Health::Hh_income_susc_mod_floor = 0.0;

double Health::health_insurance_distribution[Insurance_assignment_index::UNSET];
int Health::health_insurance_cdf_size = 0;

// static method called in main (Fred.cc)

void Health::initialize_static_variables() {
  //Setup the static variables if they aren't already initialized
  if(!Health::is_initialized) {

    Params::get_param_from_string("min_age_face_masks", &(Health::Min_age_wear_face_masks));
    Params::get_param_from_string("days_to_wear_face_masks", &(Health::Days_to_wear_face_masks));
    Params::get_param_from_string("day_start_wearing_face_masks",
				  &(Health::Day_start_wearing_face_masks));
    //Load face mask compliance map
    int N_face_mask_locations, N_face_mask_compliance;
    Params::get_param_from_string("face_mask_locations", &(N_face_mask_locations));

    Params::get_param_from_string("face_mask_compliance", &(N_face_mask_compliance));
    if (N_face_mask_compliance != N_face_mask_locations) {
      Utils::fred_abort("the number of locations and number of compliance levels don't match");
    }

    printf("N_face_mask_locations: %d and N_face_mask_compliance: %d\n", N_face_mask_locations, N_face_mask_compliance);
    char fm_str[MAX_PARAM_SIZE];
    Params::get_param((char*)"face_mask_compliance", fm_str);    
    std::vector<string> face_mask_locations_arr;
    std::vector<double> face_mask_compliance_arr;
    Params::get_param_vector_from_string((char*)"face_mask_locations", face_mask_locations_arr);
    Params::get_param_vector_from_string(fm_str,face_mask_compliance_arr);
    //Move into unordered_map
    for (int ii = 0; ii < N_face_mask_locations; ++ii) {
      printf("face_mask locations %s face_mask_compliance %lf\n", face_mask_locations_arr[ii].c_str(),
	     face_mask_compliance_arr[ii]);
      Face_mask_compliance[face_mask_locations_arr[ii]] = face_mask_compliance_arr[ii];
    }
    if (Face_mask_compliance.count("other") == 0) {
      Utils::fred_abort("there should always be a location called \"other\"");
    }

    Params::get_param_from_string("hand_washing_compliance", &(Health::Hand_washing_compliance));

    int temp_int = 0;

    if(Global::Enable_hh_income_based_susc_mod) {
      Params::get_param_from_string("hh_income_susc_mod_floor", &(Health::Hh_income_susc_mod_floor));
    }

    if(Global::Enable_Chronic_Condition) {
      Health::asthma_prob = new Age_Map("Asthma Probability");
      Health::asthma_prob->read_from_input("asthma_prob");
      Health::asthma_hospitalization_prob_mult = new Age_Map("Asthma Hospitalization Probability Mult");
      Health::asthma_hospitalization_prob_mult->read_from_input("asthma_hospitalization_prob_mult");
      Health::asthma_case_fatality_prob_mult = new Age_Map("Asthma Case Fatality Probability Mult");
      Health::asthma_case_fatality_prob_mult->read_from_input("asthma_case_fatality_prob_mult");

      Health::COPD_prob = new Age_Map("COPD Probability");
      Health::COPD_prob->read_from_input("COPD_prob");
      Health::COPD_hospitalization_prob_mult = new Age_Map("COPD Hospitalization Probability Mult");
      Health::COPD_hospitalization_prob_mult->read_from_input("COPD_hospitalization_prob_mult");
      Health::COPD_case_fatality_prob_mult = new Age_Map("COPD Case Fatality Probability Mult");
      Health::COPD_case_fatality_prob_mult->read_from_input("COPD_case_fatality_prob_mult");

      Health::chronic_renal_disease_prob = new Age_Map("Chronic Renal Disease Probability");
      Health::chronic_renal_disease_prob->read_from_input("chronic_renal_disease_prob");
      Health::chronic_renal_disease_hospitalization_prob_mult = new Age_Map("Chronic Renal Disease Hospitalization Probability Mult");
      Health::chronic_renal_disease_hospitalization_prob_mult->read_from_input("chronic_renal_disease_hospitalization_prob_mult");
      Health::chronic_renal_disease_case_fatality_prob_mult = new Age_Map("Chronic Renal Disease Case Fatality Probability Mult");
      Health::chronic_renal_disease_case_fatality_prob_mult->read_from_input("chronic_renal_disease_case_fatality_prob_mult");

      Health::diabetes_prob = new Age_Map("Diabetes Probability");
      Health::diabetes_prob->read_from_input("diabetes_prob");
      Health::diabetes_hospitalization_prob_mult = new Age_Map("Diabetes Hospitalization Probability Mult");
      Health::diabetes_hospitalization_prob_mult->read_from_input("diabetes_hospitalization_prob_mult");
      Health::diabetes_case_fatality_prob_mult = new Age_Map("Diabetes Case Fatality Probability Mult");
      Health::diabetes_case_fatality_prob_mult->read_from_input("diabetes_case_fatality_prob_mult");

      Health::heart_disease_prob = new Age_Map("Heart Disease Probability");
      Health::heart_disease_prob->read_from_input("heart_disease_prob");
      Health::heart_disease_hospitalization_prob_mult = new Age_Map("Heart Disease Hospitalization Probability Mult");
      Health::heart_disease_hospitalization_prob_mult->read_from_input("heart_disease_hospitalization_prob_mult");
      Health::heart_disease_case_fatality_prob_mult = new Age_Map("Heart Disease Case Fatality Probability Mult");
      Health::heart_disease_case_fatality_prob_mult->read_from_input("heart_disease_case_fatality_prob_mult");

      Health::hypertension_prob = new Age_Map("Hypertension Probability");
      Health::hypertension_prob->read_from_input("hypertension_prob");
      Health::hypertension_hospitalization_prob_mult = new Age_Map("Hypertension Hospitalization Probability Mult");
      Health::hypertension_hospitalization_prob_mult->read_from_input("hypertension_hospitalization_prob_mult");
      Health::hypertension_case_fatality_prob_mult = new Age_Map("Hypertension Case Fatality Probability Mult");
      Health::hypertension_case_fatality_prob_mult->read_from_input("hypertension_case_fatality_prob_mult");

      Health::hypercholestrolemia_prob = new Age_Map("Hypercholestrolemia Probability");
      Health::hypercholestrolemia_prob->read_from_input("hypercholestrolemia_prob");
      Health::hypercholestrolemia_hospitalization_prob_mult = new Age_Map("Hypercholestrolemia Hospitalization Probability Mult");
      Health::hypercholestrolemia_hospitalization_prob_mult->read_from_input("hypercholestrolemia_hospitalization_prob_mult");
      Health::hypercholestrolemia_case_fatality_prob_mult = new Age_Map("Hypercholestrolemia Case Fatality Probability Mult");
      Health::hypercholestrolemia_case_fatality_prob_mult->read_from_input("hypercholestrolemia_case_fatality_prob_mult");

      Health::pregnancy_hospitalization_prob_mult = new Age_Map("Pregnancy Hospitalization Probability Mult");
      Health::pregnancy_hospitalization_prob_mult->read_from_input("pregnancy_hospitalization_prob_mult");
      Health::pregnancy_case_fatality_prob_mult = new Age_Map("Pregnancy Case Fatality Probability Mult");
      Health::pregnancy_case_fatality_prob_mult->read_from_input("pregnancy_case_fatality_prob_mult");
    }

    if(Global::Enable_Health_Insurance) {

      Health::health_insurance_cdf_size = Params::get_param_vector((char*)"health_insurance_distribution", Health::health_insurance_distribution);

      // convert to cdf
      double stotal = 0;
      for(int i = 0; i < Health::health_insurance_cdf_size; ++i) {
        stotal += Health::health_insurance_distribution[i];
      }
      if(stotal != 100.0 && stotal != 1.0) {
        Utils::fred_abort("Bad distribution health_insurance_distribution params_str\nMust sum to 1.0 or 100.0\n");
      }
      double cumm = 0.0;
      for(int i = 0; i < Health::health_insurance_cdf_size; ++i) {
        Health::health_insurance_distribution[i] /= stotal;
        Health::health_insurance_distribution[i] += cumm;
        cumm = Health::health_insurance_distribution[i];
      }
    }

    Health::is_initialized = true;
  }
}

Insurance_assignment_index::e Health::get_health_insurance_from_distribution() {
  if(Global::Enable_Health_Insurance && Health::is_initialized) {
    int i = Random::draw_from_distribution(Health::health_insurance_cdf_size, Health::health_insurance_distribution);
    return Health::get_insurance_type_from_int(i);
  } else {
    return Insurance_assignment_index::UNSET;
  }
}

double Health::get_chronic_condition_case_fatality_prob_mult(double real_age, Chronic_condition_index::e cond_idx) {
  if(Global::Enable_Chronic_Condition && Health::is_initialized) {
    assert(cond_idx >= Chronic_condition_index::ASTHMA);
    assert(cond_idx < Chronic_condition_index::CHRONIC_MEDICAL_CONDITIONS);
    switch(cond_idx) {
    case Chronic_condition_index::ASTHMA:
      return Health::asthma_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::COPD:
      return Health::COPD_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::CHRONIC_RENAL_DISEASE:
      return Health::chronic_renal_disease_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::DIABETES:
      return Health::diabetes_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::HEART_DISEASE:
      return Health::heart_disease_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::HYPERTENSION:
      return Health::hypertension_case_fatality_prob_mult->find_value(real_age);
    case Chronic_condition_index::HYPERCHOLESTROLEMIA:
      return Health::hypercholestrolemia_case_fatality_prob_mult->find_value(real_age);
    default:
      return 1.0;
    }
  }
  return 1.0;
}

double Health::get_chronic_condition_hospitalization_prob_mult(double real_age, Chronic_condition_index::e cond_idx) {
  if(Global::Enable_Chronic_Condition && Health::is_initialized) {
    assert(cond_idx >= Chronic_condition_index::ASTHMA);
    assert(cond_idx < Chronic_condition_index::CHRONIC_MEDICAL_CONDITIONS);
    switch(cond_idx) {
    case Chronic_condition_index::ASTHMA:
      return Health::asthma_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::COPD:
      return Health::COPD_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::CHRONIC_RENAL_DISEASE:
      return Health::chronic_renal_disease_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::DIABETES:
      return Health::diabetes_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::HEART_DISEASE:
      return Health::heart_disease_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::HYPERTENSION:
      return Health::hypertension_hospitalization_prob_mult->find_value(real_age);
    case Chronic_condition_index::HYPERCHOLESTROLEMIA:
      return Health::hypercholestrolemia_hospitalization_prob_mult->find_value(real_age);
    default:
      return 1.0;
    }
  }
  return 1.0;
}

Health::Health() {
  this->myself = NULL;
  this->past_infections = NULL;
  this->alive = true;
  this->av_health = NULL;
  this->checked_for_av = NULL;
  this->vaccine_health = NULL;
  this->has_face_mask_behavior_anywhere = false;
  this->wears_face_mask_today = false;
  this->days_wearing_face_mask = 0;
  this->washes_hands = false;
  this->days_symptomatic = 0;
  this->days_hospitalization = 0;
  this->previous_infection_serotype = 0;
  this->insurance_type = Insurance_assignment_index::UNSET;
  this->infection = NULL;
  this->immunity_end_date = NULL;
  this->infectee_count = NULL;
  this->susceptibility_multp = NULL;
  this->exposure_date = NULL;
  this->infector_id = NULL;
  this->infected_in_mixing_group = NULL;
  this->health_condition = NULL;
  this->health_state.clear();
  this->total_number_of_infections = 0;
}

void Health::setup(Person* self) {
  this->myself = self;
  FRED_VERBOSE(1, "Health::setup for person %d\n", myself->get_id());
  this->alive = true;
  this->intervention_flags = intervention_flags_type();
  // infection pointers stored in statically allocated array (length of which
  // is determined by static constant Global::MAX_NUM_DISEASES)
  this->susceptible = fred::disease_bitset();
  this->infectious = fred::disease_bitset();
  this->symptomatic = fred::disease_bitset();
  this->hospitalized = fred::disease_bitset();
  this->recovered = fred::disease_bitset();
  this->immunity = fred::disease_bitset();

  // Determines if the agent is at risk
  this->at_risk = fred::disease_bitset();
  
  // Determine if the agent washes hands
  this->washes_hands = false;
  if(Health::Hand_washing_compliance > 0.0) {
    this->washes_hands = (Random::draw_random() < Health::Hand_washing_compliance);
  }
  
  // Facemasks defaults
  this->has_face_mask_behavior_anywhere = false;
  this->wears_face_mask_today = false;
  this->days_wearing_face_mask = 0;
  
  this->case_fatality = fred::disease_bitset();
  int diseases = Global::Diseases.get_number_of_diseases();
  FRED_VERBOSE(1, "Health::setup diseases %d\n", diseases);
  this->infection = new Infection* [diseases];
  this->susceptibility_multp = new double [diseases];
  this->infectee_count = new int [diseases];
  this->exposure_date = new int [diseases];
  this->infector_id = new int [diseases];
  this->infected_in_mixing_group = new Mixing_Group*[diseases];
  this->immunity_end_date = new int [diseases];
  this->past_infections = new past_infections_type [diseases];
  this->health_condition = new health_condition_t [diseases];

  for(int disease_id = 0; disease_id < diseases; ++disease_id) {
    this->recovered.reset(disease_id);
    this->susceptible.reset(disease_id);
    this->case_fatality.reset(disease_id);
    this->infection[disease_id] = NULL;
    this->susceptibility_multp[disease_id] = 1.0;
    this->infectee_count[disease_id] = 0;
    this->exposure_date[disease_id] = -1;
    this->infector_id[disease_id] = -1;
    this->infected_in_mixing_group[disease_id] = NULL;
    this->immunity_end_date[disease_id] = -1;
    this->past_infections[disease_id].clear();
    this->health_condition[disease_id].state = -1;
    this->health_condition[disease_id].last_transition_day = -1;
    this->health_condition[disease_id].next_state = -1;
    this->health_condition[disease_id].next_transition_day = -1;

    Disease* disease = Global::Diseases.get_disease(disease_id);
    if (disease->assume_susceptible()) {
      become_susceptible(disease_id);
    }

    if(disease->get_at_risk() != NULL && !disease->get_at_risk()->is_empty()) {
      double at_risk_prob = disease->get_at_risk()->find_value(myself->get_real_age());
      if(Random::draw_random() < at_risk_prob) { // Now a probability <=1.0
        declare_at_risk(disease);
      }
    }
  }

  this->days_symptomatic = 0;
  this->days_hospitalization = 0;
  this->vaccine_health = NULL;
  this->av_health = NULL;
  this->checked_for_av = NULL;
  this->previous_infection_serotype = -1;

  if(Health::nantivirals == -1) {
    Params::get_param_from_string("number_antivirals", &Health::nantivirals);
  }

  if(Global::Enable_Chronic_Condition && Health::is_initialized) {
    double prob = 0.0;
    prob = Health::asthma_prob->find_value(myself->get_real_age());
    set_is_asthmatic((Random::draw_random() < prob));

    prob = Health::COPD_prob->find_value(myself->get_real_age());
    set_has_COPD((Random::draw_random() < prob));

    prob = Health::chronic_renal_disease_prob->find_value(myself->get_real_age());
    set_has_chronic_renal_disease((Random::draw_random() < prob));

    prob = Health::diabetes_prob->find_value(myself->get_real_age());
    set_is_diabetic((Random::draw_random() < prob));

    prob = Health::heart_disease_prob->find_value(myself->get_real_age());
    set_has_heart_disease((Random::draw_random() < prob));

    prob = Health::hypertension_prob->find_value(myself->get_real_age());
    set_has_hypertension((Random::draw_random() < prob));

    prob = Health::hypercholestrolemia_prob->find_value(myself->get_real_age());
    set_has_hypercholestrolemia((Random::draw_random() < prob));
  }
}

Health::~Health() {
  for(size_t i = 0; i < Global::Diseases.get_number_of_diseases(); ++i) {
    if(this->infection[i] != NULL) {
      delete this->infection[i];
    }
  }

  if(this->vaccine_health) {
    for(unsigned int i = 0; i < this->vaccine_health->size(); ++i) {
      delete (*this->vaccine_health)[i];
    }
    this->vaccine_health->clear();
    delete this->vaccine_health;
  }

  if(this->av_health) {
    for(unsigned int i = 0; i < this->av_health->size(); ++i) {
      delete (*this->av_health)[i];
    }
    this->av_health->clear();
    delete this->av_health;
  }

  if(this->checked_for_av) {
    delete this->checked_for_av;
  }
}

void Health::become_susceptible(int disease_id) {
  if(this->susceptible.test(disease_id)) {
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			     "HEALTH CHART: %s person %d is already SUSCEPTIBLE for disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
    return;
  }  
  if(this->infection[disease_id] != NULL){
    printf("HEALTH CHART: %s person %d has an active infection for disease %d\n",
	   Date::get_date_string().c_str(),
	   myself->get_id(), disease_id);
  }
  assert(this->infection[disease_id] == NULL);
  this->susceptibility_multp[disease_id] = 1.0;
  this->susceptible.set(disease_id);

  assert(is_susceptible(disease_id));
  this->recovered.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is SUSCEPTIBLE for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::become_susceptible_by_natural_waning(int disease_id) {
  if(this->susceptible.test(disease_id)) {
    this->immunity.reset(disease_id);
    return;
  }  
  if(this->infection[disease_id] == NULL) {
    // not already infected
    this->susceptibility_multp[disease_id] = 1.0;

    this->susceptible.set(disease_id);
    this->immunity.reset(disease_id);
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d is SUSCEPTIBLE for disease %d by natural waning\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
  } else {
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d had no vaccine waning because was already infected with disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
  }
}

void Health::become_susceptible_by_vaccine_waning(int disease_id) {
  if(this->susceptible.test(disease_id)) {
    return;
  }
  if(this->infection[disease_id] == NULL) {
    // not already infected
    this->susceptibility_multp[disease_id] = 1.0;
    this->susceptible.set(disease_id);
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d is SUSCEPTIBLE for disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
  } else {
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d had no vaccine waning because was already infected with disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
  }
}

void Health::become_susceptible_to_symptoms_by_vaccine_waning(int disease_id) {
  this->immunity_to_symptoms.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			   "HEALTH CHART: %s person %d is SUSCEPTIBLE TO SYMPTOMS for disease %d\n",
			   Date::get_date_string().c_str(),
			   myself->get_id(), disease_id);
}

void Health::become_susceptible_to_hospitalization_by_vaccine_waning(int disease_id) {
  this->immunity_to_hospitalization.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			   "HEALTH CHART: %s person %d is SUSCEPTIBLE TO HOSPITALIZATION for disease %d\n",
			   Date::get_date_string().c_str(),
			   myself->get_id(), disease_id);
}

void Health::become_exposed(int disease_id, Person* infector, Mixing_Group* mixing_group, int day) {

  FRED_VERBOSE(0, "become_exposed: person %d is exposed to disease %d day %d\n",
		            myself->get_id(), disease_id, day);
  // if(myself->get_id() == 7147928){
  //   printf("become_exposed: person %d is exposed to disease %d day %d\n",
  // 	   myself->get_id(), disease_id, day);
  // }
  if(this->infection[disease_id] != NULL) {
    Utils::fred_abort("DOUBLE EXPOSURE: person %d dis_id %d day %d\n", myself->get_id(), disease_id, day);
  }

  if(Global::Verbose > 0) {
    if(mixing_group == NULL) {
      FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			      "HEALTH CHART: %s person %d is an IMPORTED EXPOSURE to disease %d\n",
			      Date::get_date_string().c_str(),
			      myself->get_id(), disease_id);
    } else {
      FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			      "HEALTH CHART: %s person %d is EXPOSED to disease %d\n",
			      Date::get_date_string().c_str(),
			      myself->get_id(), disease_id);
    }
  }
  this->total_number_of_infections++;
  this->infectious.reset(disease_id);
  this->symptomatic.reset(disease_id);
  this->hospitalized.reset(disease_id);
  Disease* disease = Global::Diseases.get_disease(disease_id);
  this->infection[disease_id] = Infection::get_new_infection(disease, infector, myself, mixing_group, day);
  FRED_VERBOSE(1, "setup infection: person %d dis_id %d day %d\n", myself->get_id(), disease_id, day);
  this->infection[disease_id]->setup();
  this->infection[disease_id]->report_infection(day);
  become_unsusceptible(disease);
  this->immunity_end_date[disease_id] = -1;
  if(myself->get_household() != NULL) {
    myself->get_household()->set_exposed(disease_id);
    myself->set_exposed_household(myself->get_household()->get_index());
  }
  if(infector != NULL) {
    this->infector_id[disease_id] = infector->get_id();
  }
  this->exposure_date[disease_id] = day;
  this->infected_in_mixing_group[disease_id] = mixing_group;

  /*
  if(myself->get_id() == 3936721){
    printf("PERSON ID: infectious %d-%d %d symptoms %d to %d hospitalization %d to %d\n", myself->get_id(),
	   this->get_infectious_start_date(disease_id), this->get_infectious_end_date(disease_id),
	   this->get_symptoms_start_date(disease_id), this->get_symptoms_end_date(disease_id),
	   this->get_hospitalization_start_date(disease_id), this->get_hospitalization_end_date(disease_id));
    abort();
    }
  */
  if(Global::Enable_Transmission_Network) {
    FRED_VERBOSE(1, "Joining transmission network: %d\n", myself->get_id());
    myself->join_network(Global::Transmission_Network);
  }
  if(Global::Enable_Disease_Cross_Protection == true && Global::Diseases.get_number_of_diseases() > 1){
    for(int dis_i = 0; dis_i < Global::Diseases.get_number_of_diseases(); ++dis_i){
      if(dis_i == disease_id){
	continue;
      }
      this->infectious.reset(dis_i);
      this->symptomatic.reset(dis_i);
      this->hospitalized.reset(dis_i);
      this->immunity_end_date[dis_i] = -1;
      this->immunity.reset(dis_i);
      become_unsusceptible(Global::Diseases.get_disease(dis_i));      
    }
  }
  if(Global::Enable_Vector_Transmission && Global::Diseases.get_number_of_diseases() > 1) {
    // special check for multi-serotype dengue:
    if(this->previous_infection_serotype == -1) {
      // remember this infection's serotype
      this->previous_infection_serotype = disease_id;
      // after the first infection, become immune to other two serotypes.
      for(int sero = 0; sero < Global::Diseases.get_number_of_diseases(); ++sero) {
        // if (sero == previous_infection_serotype) continue;
        if(sero == disease_id) {
          continue;
        }
        FRED_STATUS(1, "DENGUE: person %d now immune to serotype %d\n",
		    myself->get_id(), sero);
        become_unsusceptible(Global::Diseases.get_disease(sero));
      }
    } else {
      // after the second infection, become immune to other two serotypes.
      // This isn't making everyone immune, just unsusceptible
      for(int sero = 0; sero < Global::Diseases.get_number_of_diseases(); ++sero) {
        if(sero == this->previous_infection_serotype) {
          continue;
        }
        if(sero == disease_id) {
          continue;
        }
        FRED_STATUS(1, "DENGUE: person %d now immune to serotype %d\n",
		    myself->get_id(), sero);
        become_unsusceptible(Global::Diseases.get_disease(sero));
      }
    }
  }
}

void Health::become_unsusceptible(int disease_id) {
  this->susceptible.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is UNSUSCEPTIBLE for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::become_unsusceptible(Disease* disease) {
  int disease_id = disease->get_id();
  become_unsusceptible(disease_id);
}

void Health::become_infectious(Disease* disease) {
  int disease_id = disease->get_id();
  assert(this->infection[disease_id] != NULL);
  this->infectious.set(disease_id);
  int household_index = myself->get_exposed_household_index();
  Household* h = Global::Places.get_household_ptr(household_index);
  assert(h != NULL);
  h->set_human_infectious(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is INFECTIOUS for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);

  /*
    If facemasks enabled, then decide if will wear facemasks 
  */

  // Determine if the agent will wear a face mask
  if(Global::Enable_Face_Mask_Usage == true && myself->get_age() >= Health::Min_age_wear_face_masks){
    if(Global::Enable_Face_Mask_Timeseries_File == true){
      for (auto it = Face_mask_compliance.begin(); it != Face_mask_compliance.end(); ++it) {
	double tmp_compliance = Global::Places.get_face_mask_compliance_today(it->first);
	//printf("Health.cc::BECOME_INFECTIOUS. Compliance %lf location %s\n", tmp_compliance, it->first.c_str());
	if(tmp_compliance > 0.0 && Random::draw_random() < tmp_compliance){
	  this->has_face_mask_behavior[it->first] = true;
	  this->has_face_mask_behavior_anywhere = true;
	}else{
	  this->has_face_mask_behavior[it->first] = false;
	}
      }
    }else{
      for (auto it = Face_mask_compliance.begin(); it != Face_mask_compliance.end(); ++it) {  
	if((it->second > 0.0) && Random::draw_random() < it->second) {
	  this->has_face_mask_behavior[it->first] = true;
	  this->has_face_mask_behavior_anywhere = true;
	} else {
	  this->has_face_mask_behavior[it->first] = false;
	}
      }
    }
  }
}

void Health::become_noninfectious(Disease* disease) {
  int disease_id = disease->get_id();
  assert(this->infection[disease_id] != NULL);
  this->infectious.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is NONINFECTIOUS for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::become_symptomatic(Disease* disease) {
  int disease_id = disease->get_id();
  if(this->infection[disease_id] == NULL) {
    FRED_STATUS(1, "Help: becoming symptomatic with no infection: person %d, disease_id %d\n", myself->get_id(), disease_id);
  }
  assert(this->infection[disease_id] != NULL);
  if(this->symptomatic.test(disease_id)) {
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d is ALREADY SYMPTOMATIC for disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
    return;
  }
  this->symptomatic.set(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is SYMPTOMATIC for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::become_hospitalized(Disease* disease) {
  int disease_id = disease->get_id();
  if(this->infection[disease_id] == NULL) {
    printf("Help: becoming hospitalized with no infection: person %d, disease_id %d\n", myself->get_id(), disease_id);
    FRED_STATUS(1, "Help: becoming hospitalized with no infection: person %d, disease_id %d\n", myself->get_id(), disease_id);
  }
  assert(this->infection[disease_id] != NULL);
  if(this->hospitalized.test(disease_id)) {
    FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			    "HEALTH CHART: %s person %d is ALREADY HOSPITALIZED for disease %d\n",
			    Date::get_date_string().c_str(),
			    myself->get_id(), disease_id);
    return;
  }
  this->hospitalized.set(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is HOSPITALIZED for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::resolve_symptoms(Disease* disease) {
  int disease_id = disease->get_id();
  // assert(this->infection[disease_id] != NULL);
  if(this->symptomatic.test(disease_id)) {
    this->symptomatic.reset(disease_id);
  }
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d RESOLVES SYMPTOMS for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}


void Health::resolve_hospitalization(Disease* disease) {
  int disease_id = disease->get_id();
  // assert(this->infection[disease_id] != NULL);
  if(this->hospitalized.test(disease_id)) {
    this->hospitalized.reset(disease_id);
  }
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d RESOLVES HOSPITALIZATION for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}


void Health::recover(Disease* disease, int day) {
  int disease_id = disease->get_id();
  // assert(this->infection[disease_id] != NULL);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			   "HEALTH CHART: %s person %d is RECOVERED from disease %d\n",
			   Date::get_date_string().c_str(),
			   myself->get_id(), disease_id);
  this->recovered.set(disease_id);
  int household_index = myself->get_exposed_household_index();
  Household* h = Global::Places.get_household_ptr(household_index);
  h->set_recovered(disease_id);
  h->reset_human_infectious();
  myself->reset_neighborhood();

  this->immunity_end_date[disease_id] = this->infection[disease_id]->get_immunity_end_date();
  if (this->immunity_end_date[disease_id] > -1) {
    this->immunity_end_date[disease_id] += day;
  }
  become_removed(disease_id, day);
}

void Health::become_removed(int disease_id, int day) {
  terminate_infection(disease_id, day);
  this->susceptible.reset(disease_id);
  this->infectious.reset(disease_id);
  this->symptomatic.reset(disease_id);
  this->hospitalized.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			   "HEALTH CHART: %s person %d is REMOVED for disease %d\n",
			   Date::get_date_string().c_str(),
			   myself->get_id(), disease_id);
}

void Health::become_immune(Disease* disease) {
  int disease_id = disease->get_id();
  disease->become_immune(myself, this->susceptible.test(disease_id),
			 this->infectious.test(disease_id), this->symptomatic.test(disease_id), this->hospitalized.test(disease_id));
  this->immunity.set(disease_id);
  this->susceptible.reset(disease_id);
  this->infectious.reset(disease_id);
  this->symptomatic.reset(disease_id);
  this->hospitalized.reset(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is IMMUNE for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}

void Health::become_immune_to_symptoms(Disease* disease) {
  int disease_id = disease->get_id();
  disease->become_immune_to_symptoms(myself, this->susceptible.test(disease_id),
				     this->infectious.test(disease_id), this->symptomatic.test(disease_id),this->hospitalized.test(disease_id));
  this->immunity_to_symptoms.set(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is IMMUNE TO SYMPTOMS for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}


void Health::become_immune_to_hospitalization(Disease* disease) {
  int disease_id = disease->get_id();
  disease->become_immune_to_hospitalization(myself, this->susceptible.test(disease_id),
					    this->infectious.test(disease_id), this->symptomatic.test(disease_id), this->hospitalized.test(disease_id));
  this->immunity_to_hospitalization.set(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			  "HEALTH CHART: %s person %d is IMMUNE TO HOSPITALIZATION for disease %d\n",
			  Date::get_date_string().c_str(),
			  myself->get_id(), disease_id);
}


void Health::become_case_fatality(int disease_id, int day) {
  FRED_VERBOSE(0, "DISEASE %d is FATAL: day %d person %d\n", disease_id, day, myself->get_id());
  this->case_fatality.set(disease_id);
  FRED_CONDITIONAL_VERBOSE(0, Global::Enable_Health_Charts,
			   "HEALTH CHART: %s person %d is CASE_FATALITY for disease %d\n",
			   Date::get_date_string().c_str(),
			   myself->get_id(), disease_id);

  /* EDIT: We have to TERMINATE a fatality case before become_removed() because 
     become removed cleans all the disease information...
   */
  
  //Have to figure out a way to remove all events from all diseases when people die
  //Right now if a person dies from one disease, the other diseases don't update
  Global::Diseases.get_disease(disease_id)->terminate_person(myself, day);
  
  
  become_removed(disease_id, day);

  // update household counts
  Mixing_Group* hh = myself->get_household();
  if(hh == NULL) {
    if(Global::Enable_Hospitals && myself->is_hospitalized() && myself->get_permanent_household() != NULL) {
      hh = myself->get_permanent_household();
    }
  }
  if(hh != NULL) {
    hh->increment_case_fatalities(day, disease_id);
  }


  // queue removal from population
  Global::Pop.prepare_to_die(day, myself);
}

void Health::update_infection(int day, int disease_id) {
  
  if(this->has_face_mask_behavior_anywhere) {
    update_face_mask_decision(day, disease_id);
  }
  
  if(this->infection[disease_id] == NULL) {
    return;
  }
  
  FRED_VERBOSE(1, "update_infection %d on day %d person %d\n", disease_id, day, myself->get_id());
  this->infection[disease_id]->update(day);
  
  // update days_symptomatic if needed
  if(this->is_symptomatic(disease_id)) {
    int days_symp_so_far = (day - this->get_symptoms_start_date(disease_id));
    if(days_symp_so_far > this->days_symptomatic) {
      this->days_symptomatic = days_symp_so_far;
    }
  }

  // update days_symptomatic if needed
  if(this->is_hospitalized(disease_id)) {
    int days_hosp_so_far = (day - this->get_hospitalization_start_date(disease_id));
    if(days_hosp_so_far > this->days_hospitalization) {
      this->days_hospitalization = days_hosp_so_far;
    }
  }

  
  // case_fatality?
  if(this->infection[disease_id]->is_fatal(day)) {
    become_case_fatality(disease_id, day);
  }
  
  FRED_VERBOSE(1,"update_infection %d FINISHED on day %d person %d\n",
	       disease_id, day, myself->get_id());

} // end Health::update_infection //


void Health::update_face_mask_decision(int day, int disease_id) {
  // printf("update_face_mask_decision entered on day %d for person %d\n", day, myself->get_id());
  Disease* disease = Global::Diseases.get_disease(disease_id);

  /*
    IF TIMESTEP IS NOT BEING READ, then use basic model, 
    IF TIMESTEP BEING READ, then go to Global::Places and get the daily probability
  */
  
  // should we start use face mask?
  if((!disease->get_face_mask_symptomatic_only() || this->is_symptomatic(day))
     && this->days_wearing_face_mask == 0
     && this->Day_start_wearing_face_masks <= day) {
    FRED_VERBOSE(1, "FACEMASK: person %d starts wearing face mask on day %d\n", myself->get_id(), day);
    this->start_wearing_face_mask();
  }

  // should we stop using face mask?
  if(this->is_wearing_face_mask()) {
    if ((!disease->get_face_mask_symptomatic_only() || this->is_symptomatic(day))
	&& this->days_wearing_face_mask < Health::Days_to_wear_face_masks) {
      this->days_wearing_face_mask++;
    } else {
      FRED_VERBOSE(1, "FACEMASK: person %d stops wearing face mask on day %d\n", myself->get_id(), day);
      this->stop_wearing_face_mask();
    }
  }
}

void Health::update_vaccine_interventions(int day){
  if(!(this->alive)) {
    return;
  }
  if(this->intervention_flags[Intervention_flag::TAKES_VACCINE]) {
    int size = (int)(this->vaccine_health->size());
    for(int i = 0; i < size; ++i) {
      (*this->vaccine_health)[i]->update(day, myself->get_real_age());
    }
  }  
}

int Health::is_vaccine_effective_any() const{
  if(this->vaccine_health){
    return (*this->vaccine_health)[0]->is_effective_any();
  }else{
    return -1;
  }
}

int Health::get_vaccination_any_effective_day() const{
  if(this->vaccine_health){
    return (*this->vaccine_health)[0]->get_vaccination_any_effective_day();
    }else{
    return -1;
  }
}

int Health::get_vaccination_immunity_loss_day() const {
  if(this->vaccine_health){
    return (*this->vaccine_health)[0]->get_vaccine_immunity_loss_day();
  }else{
    return -1;
  }
}

void Health::update_interventions(int day) {
  // if deceased, health status should have been cleared during population
  // update (by calling Person->die(), then Health->die(), which will reset (bool) alive
  if(!(this->alive)) {
    return;
  }
  //if(this->intervention_flags.any()) {
  // update antiviral status
  if(this->intervention_flags[Intervention_flag::TAKES_AV]) {
    for(av_health_itr i = this->av_health->begin(); i != this->av_health->end(); ++i) {
      (*i)->update(day);
    }
  }
} // end Health::update_interventions

void Health::declare_at_risk(Disease* disease) {
  int disease_id = disease->get_id();
  this->at_risk.set(disease_id);
}

void Health::advance_seed_infection(int disease_id, int days_to_advance) {
  assert(this->infection[disease_id] != NULL);
  this->infection[disease_id]->advance_seed_infection(days_to_advance);
}

int Health::get_exposure_date(int disease_id) const {
  return this->exposure_date[disease_id];
}

int Health::get_infectious_start_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_infectious_start_date();
  }
}

int Health::get_infectious_end_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_infectious_end_date();
  }
}

int Health::get_symptoms_start_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_symptoms_start_date();
  }
}

int Health::get_symptoms_end_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_symptoms_end_date();
  }
}

int Health::get_hospitalization_start_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_hospitalization_start_date();
  }
}

int Health::get_hospitalization_end_date(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return -1;
  } else {
    return this->infection[disease_id]->get_hospitalization_end_date();
  }
}

int Health::get_immunity_end_date(int disease_id) const {
  return this->immunity_end_date[disease_id];
}

bool Health::is_recovered(int disease_id) {
  return this->recovered.test(disease_id);
}


int Health::get_infector_id(int disease_id) const {
  return infector_id[disease_id];
}

Person* Health::get_infector(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    return NULL;
  } else {
    return this->infection[disease_id]->get_infector();
  }
}

Mixing_Group* Health::get_infected_mixing_group(int disease_id) const {
  return this->infected_in_mixing_group[disease_id];
}

int Health::get_infected_mixing_group_id(int disease_id) const {
  Mixing_Group* mixing_group = get_infected_mixing_group(disease_id);
  if(mixing_group == NULL) {
    return -1;
  } else {
    return mixing_group->get_id();
  }
}

char Health::get_infected_mixing_group_type(int disease_id) const {
  Mixing_Group* mixing_group = get_infected_mixing_group(disease_id);
  if(mixing_group == NULL) {
    return 'X';
  } else {
    return mixing_group->get_type();
  }
}

char dummy_label[8];
char* Health::get_infected_mixing_group_label(int disease_id) const {
  if(this->infection[disease_id] == NULL) {
    strcpy(dummy_label, "-");
    return dummy_label;
  }
  Mixing_Group* mixing_group = get_infected_mixing_group(disease_id);
  if(mixing_group == NULL) {
    strcpy(dummy_label, "X");
    return dummy_label;
  } else {
    return mixing_group->get_label();
  }
}

int Health::get_infectees(int disease_id) const {
  return this->infectee_count[disease_id];
}

double Health::get_susceptibility(int disease_id) const {
  double suscep_multp = this->susceptibility_multp[disease_id];

  if(this->infection[disease_id] == NULL) {
    return suscep_multp;
  } else {
    return this->infection[disease_id]->get_susceptibility() * suscep_multp;
  }
}

double Health::get_infectivity(int disease_id, int day) const {
  if(this->infection[disease_id] == NULL) {
    return 0.0;
  } else {
    return this->infection[disease_id]->get_infectivity(day);
  }
}

double Health::get_shedding(int disease_id, int day) const {
  if(this->infection[disease_id] == NULL) {
    return 0.0;
  } else {
    return this->infection[disease_id]->get_shedding(day);
  }
}

double Health::get_symptoms(int disease_id, int day) const {

  if(this->infection[disease_id] == NULL) {
    return 0.0;
  } else {
    return this->infection[disease_id]->get_symptoms(day);
  }
}

//Modify Operators

double Health::get_susceptibility_modifier_due_to_person_age(int disease_id, int int_age) {
  Disease* disease = Global::Diseases.get_disease(disease_id);
  if(Global::Enable_Age_Specific_Susceptibility){
    return (disease->get_age_specific_susceptibility_efficacy(int_age));
  }
  return 1.0;
}

double Health::get_transmission_modifier_due_to_hygiene(int disease_id, Place* place) {
  Disease* disease = Global::Diseases.get_disease(disease_id);
  bool face_mask_here = false;
  if (disease->is_face_mask_usage_enabled()) {
    if (this->is_wearing_face_mask()) {
      string loc_type(place->get_place_type());
      string loc_subtype(place->get_place_subtype());
      string loc_hhtype(place->get_household_type());
      if (this->has_face_mask_behavior.count(loc_hhtype) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_hhtype];
      } else if (this->has_face_mask_behavior.count(loc_subtype) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_subtype];
      } else if (this->has_face_mask_behavior.count(loc_type) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_type];
      } else {
	face_mask_here = this->has_face_mask_behavior["other"];
      }
    }
  }
  if(face_mask_here && disease->is_hand_washing_enabled() && this->is_washing_hands()) {
    return (1.0 - disease->get_face_mask_plus_hand_washing_transmission_efficacy());
  }
  if(face_mask_here) {
    if(!disease->get_face_mask_odds_ratio_method()) {
      return (1.0 - disease->get_face_mask_transmission_efficacy());
    }
  }
  if(disease->is_hand_washing_enabled() && this->is_washing_hands()) {
    return (1.0 - disease->get_hand_washing_transmission_efficacy());
  }
  return 1.0;
}

double Health::get_susceptibility_modifier_due_to_hygiene(int disease_id) {
  Disease* disease = Global::Diseases.get_disease(disease_id);
  if(disease->is_hand_washing_enabled() && this->is_washing_hands()) {
    return (1.0 - disease->get_hand_washing_susceptibility_efficacy());
  }
  return 1.0;
}

double Health::get_infection_modifier_face_masks_odds_ratio(int disease_id, double infection_prob,
							    Place* place) {
  Disease* disease = Global::Diseases.get_disease(disease_id);
  bool face_mask_here = false;
  if (disease->is_face_mask_usage_enabled()) {
    if (this->is_wearing_face_mask()) {
      string loc_type(place->get_place_type());
      string loc_subtype(place->get_place_subtype());
      string loc_hhtype(place->get_household_type());
      if (this->has_face_mask_behavior.count(loc_hhtype) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_hhtype];
      } else if (this->has_face_mask_behavior.count(loc_subtype) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_subtype];
      } else if (this->has_face_mask_behavior.count(loc_type) == 1) {
	face_mask_here = this->has_face_mask_behavior[loc_type];
      } else {
	face_mask_here = this->has_face_mask_behavior["other"];
      }
    }
  }
  if(face_mask_here) {
    return (1/((1-infection_prob)/disease->get_face_mask_transmission_efficacy() + infection_prob));
  }
  return 1.0;
}


double Health::get_susceptibility_modifier_due_to_household_income(int hh_income) {

  if(Global::Enable_hh_income_based_susc_mod) {
    if(hh_income >= Household::get_min_hh_income_90_pct()) {
      return Health::Hh_income_susc_mod_floor;
    } else {
      double rise = 1.0 - Health::Hh_income_susc_mod_floor;
      double run = static_cast<double>(Household::get_min_hh_income() - Household::get_min_hh_income_90_pct());
      double m = rise / run;

      // Equation of line is y - y1 = m(x - x1)
      // y = m*x - m*x1 + y1
      double x = static_cast<double>(hh_income);
      return m * x - m * Household::get_min_hh_income() + 1.0;
    }
  } else {
    return 1.0;
  }
}

void Health::modify_susceptibility(int disease_id, double multp) {
  this->susceptibility_multp[disease_id] *= multp;
}

void Health::modify_infectivity(int disease_id, double multp) {
  if(this->infection[disease_id] != NULL) {
    this->infection[disease_id]->modify_infectivity(multp);
  }
}

void Health::modify_infectious_period(int disease_id, double multp, int cur_day) {
  if(this->infection[disease_id] != NULL) {
    this->infection[disease_id]->modify_infectious_period(multp, cur_day);
  }
}

void Health::modify_asymptomatic_period(int disease_id, double multp, int cur_day) {
  if(this->infection[disease_id] != NULL) {
    this->infection[disease_id]->modify_asymptomatic_period(multp, cur_day);
  }
}

void Health::modify_symptomatic_period(int disease_id, double multp, int cur_day) {
  if(this->infection[disease_id] != NULL) {
    this->infection[disease_id]->modify_symptomatic_period(multp, cur_day);
  }
}

void Health::modify_develops_symptoms(int disease_id, bool symptoms, int cur_day) {
  if(this->infection[disease_id] != NULL
     && ((this->infection[disease_id]->is_infectious(cur_day)
	        && !this->infection[disease_id]->is_symptomatic(cur_day))
	        || !this->infection[disease_id]->is_infectious(cur_day))) {

    this->infection[disease_id]->modify_develops_symptoms(symptoms, cur_day);
    this->symptomatic.set(disease_id);
  }
}
int Health::get_current_vaccine_dose(int i){
  if(this->vaccine_health) {
    return (*this->vaccine_health)[i]->get_current_dose();
  } else {
    return -1;
  }
}

int Health::get_days_to_next_dose(int i){
   if(this->vaccine_health) {
    return (*this->vaccine_health)[i]->get_days_to_next_dose();
  } else {
    return -1;
  }
}

int Health::get_next_dose_mix_match(int i){
   if(this->vaccine_health) {
    return (*this->vaccine_health)[i]->get_next_dose_mix_match();
  } else {
    return -1;
  }
}

//Medication operators
void Health::take_vaccine(Vaccine* vaccine, int day, Vaccine_Manager* vm) {
  // Compliance will be somewhere else
  double real_age = myself->get_real_age();
  // Is this our first dose?
  Vaccine_Health * vaccine_health_for_dose = NULL;

  if(this->vaccine_health == NULL) {
    this->vaccine_health = new vaccine_health_type();
  }

  for(unsigned int ivh = 0; ivh < this->vaccine_health->size(); ++ivh) {
    if((*this->vaccine_health)[ivh]->get_vaccine() == vaccine) {
      vaccine_health_for_dose = (*this->vaccine_health)[ivh];
    }
  }
  
  if(vaccine_health_for_dose == NULL) { // This is our first dose of this vaccine
    this->vaccine_health->push_back(new Vaccine_Health(day, vaccine, real_age, myself, vm));
    this->intervention_flags[Intervention_flag::TAKES_VACCINE] = true;
  } else { // Already have a dose, need to take the next dose
    vaccine_health_for_dose->update_for_next_dose(day, real_age);
  }

  if(Global::VaccineTracefp != NULL) {
    fprintf(Global::VaccineTracefp, " id %7d vaccid %3d vaccine health %d", myself->get_id(),
	    (*this->vaccine_health)[this->vaccine_health->size() - 1]->get_vaccine()->get_ID(), this->vaccine_health->size());
    (*this->vaccine_health)[this->vaccine_health->size() - 1]->printTrace();
    fprintf(Global::VaccineTracefp, "\n");
  }
  
  return;
}

void Health::take(Antiviral* av, int day) {
  if(this->checked_for_av == NULL) {
    this->checked_for_av = new checked_for_av_type();
    this->checked_for_av->assign(nantivirals, false);
  }
  if(this->av_health == NULL) {
    this->av_health = new av_health_type();
  }
  this->av_health->push_back(new AV_Health(day, av, this));
  this->intervention_flags[Intervention_flag::TAKES_AV] = true;
  return;
}

bool Health::is_on_av_for_disease(int day, int d) const {
  for(unsigned int iav = 0; iav < this->av_health->size(); ++iav) {
    if((*this->av_health)[iav]->get_disease() == d
       && (*this->av_health)[iav]->is_on_av(day)) {
      return true;
    }
  }
  return false;
}

int Health::get_av_start_day(int i) const {
  assert(this->av_health != NULL);
  return (*this->av_health)[i]->get_av_start_day();
}

void Health::infect(Person* infectee, int disease_id, Mixing_Group* mixing_group, int day) {
  infectee->become_exposed(disease_id, myself, mixing_group, day);
  
#pragma omp atomic
  ++(this->infectee_count[disease_id]);
  
  int exp_day = this->get_exposure_date(disease_id);
  assert(0 <= exp_day);
  Disease* disease = Global::Diseases.get_disease(disease_id);
  disease->increment_cohort_infectee_count(exp_day);

  FRED_STATUS(1, "person %d infected person %d infectees = %d\n",
	      myself->get_id(), infectee->get_id(), infectee_count[disease_id]);

  if(Global::Enable_Transmission_Network) {
    FRED_VERBOSE(1, "Creating link in transmission network: %d -> %d\n", myself->get_id(), infectee->get_id());
    myself->create_network_link_to(infectee, Global::Transmission_Network);
  }
}

void Health::update_mixing_group_counts(int day, int disease_id, Mixing_Group* mixing_group) {
  // this is only called for people with an active infection
  assert(is_infected(disease_id));

  // mixing group must exist to update
  if(mixing_group == NULL) {
    return;
  }

  // update infection counters
  if(is_newly_infected(day, disease_id)) {
    mixing_group->increment_new_infections(day, disease_id);
  }
  mixing_group->increment_current_infections(day, disease_id);
  
  // update symptomatic infection counters
  if(is_symptomatic(disease_id)) {
    if(is_newly_symptomatic(day, disease_id)) {
      mixing_group->increment_new_symptomatic_infections(day, disease_id);
    }
    mixing_group->increment_current_symptomatic_infections(day, disease_id);
  }

  // update preinfectious counters
  if (get_infectious_start_date(disease_id) > day) {
    mixing_group->increment_current_preinfectious_people(day, disease_id);
  }
}

void Health::terminate_infection(int disease_id, int day) {
  if(this->health_condition[disease_id].state > -1) {
    Global::Diseases.get_disease(disease_id)->terminate_person(myself, day);
  }
  if(this->infection[disease_id] != NULL) {
    // delete the infection object
    delete this->infection[disease_id];
    this->infection[disease_id] = NULL;
  }
}

void Health::terminate(int day) {
  for(int disease_id = 0; disease_id < Global::Diseases.get_number_of_diseases(); ++disease_id) {
    if(this->infection[disease_id] != NULL) {
      become_removed(disease_id, day);
    }
    if(this->health_condition[disease_id].state == 0) {
      Global::Diseases.get_disease(disease_id)->terminate_person(myself, day);
    }
  }
  this->alive = false;
}

void Health::update_health_conditions(int day) {
  for(int disease_id = 0; disease_id < Global::Diseases.get_number_of_diseases(); ++disease_id) {
    if(this->health_condition[disease_id].state > -1) {
      Global::Diseases.get_disease(disease_id)->get_epidemic()->transition_person(this->myself, day, this->health_condition[disease_id].state);
    }    
  }  
}

int Health::get_vaccinated_id() const {
  if(this->vaccine_health){
    return (*this->vaccine_health)[0]->get_vaccine()->get_ID();
  }else{
    return -1;
  }
}
