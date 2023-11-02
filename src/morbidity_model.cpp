// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <Rmath.h>
#include <RcppNumerical.h>
#include <cmath>

using namespace Numer;
using namespace Rcpp;
using namespace std;

// define global hard-wired parameters
const int na = 51, amax = 51;
const double da = na/amax, aged=1/da;
const int nlg = 3; // number of liver groups (AB, C, DEF)
const double shape = 1.11, scale = 20.69;

// define parameters that are assigned values at run time
double tau1, tau2, eta1, eta2;
int nsteps;

//storage vectors for output from time & age specific lists from transmission model
NumericVector cur_wb_female(na), cur_wb_male(na);  //current worm burden 
NumericVector cum_wb_female(na), cum_wb_male(na);  //cumulative worm burden 

// storage vectors for output 
NumericVector ab_age_female(na),  ab_age_male(na);
NumericVector c_age_female(na),   c_age_male(na);
NumericVector def_age_female(na), def_age_male(na);
double ab_pop_female, ab_pop_male, c_pop_female, c_pop_male, def_pop_female, def_pop_male; //population mean proportions of different liver scores

NumericVector AB_C_female(na), C_AB_female(na), AB_C_ratio_female(na), 
C_DEF_female(na), DEF_C_female(na), C_DEF_ratio_female(na);
double AB_C_ratio_pop_female, C_DEF_ratio_pop_female;

// define global storage vectors and scales
NumericVector age(na);                                                                      // age mid points
NumericVector demographic_weights(na);

// define matrix of current derivatives (simplify to start, non-sex-structured)
NumericMatrix dymat_female(na, nlg), dymat_male(na, nlg);  //3 columns for AB, C, & DEF liver fibrosis scores

// define matrices for RK4 integration
NumericMatrix k1_female(na, nlg), k1_male(na, nlg), k2_female(na, nlg), k2_male(na, nlg), k3_female(na, nlg), k3_male(na, nlg), k4_female(na, nlg), k4_male(na, nlg);

// define state variable matrices
NumericMatrix female_states(na, nlg), male_states(na, nlg), female_updated_states(na, nlg), male_updated_states(na, nlg);
NumericMatrix female_tmp_states(na, nlg), male_tmp_states(na, nlg), female_init_states(na, nlg), male_init_states(na, nlg);

// function to calculate age midpoints
void AgeMidpoints() {
  for (int i = 0; i < na; i++) 
    age[i] = i*aged + 0.5*aged;
}

// function to calculate demographic weights
void DemographicWeights() {
  for(int i = 0; i < na; i++) {
    double ageplus = age[i] + 0.5*aged;
    double ageminus = age[i] - 0.5*aged;
    demographic_weights[i] = (R::pweibull(ageplus, shape, scale,1,0) - 
      R::pweibull(ageminus, shape, scale,1,0)) / R::pweibull(amax,  shape,  scale,1,0);
  }
}

// function to initialise state variables
void InitialiseStates() {
  for(int i = 0; i < na; i++) {
    female_states(i,0) = 1; //all initially in AB state (healthy livers)
    female_states(i,1) = 0;
    female_states(i,2) = 0;
    male_states(i,0) = 1; //all initially in AB state (healthy livers)
    male_states(i,1) = 0;
    male_states(i,2) = 0;
  }
}

// function to initialize liver scores
void InitialiseLiverScores() {
  for (int i=0; i<na; i++) {
    ab_age_female[i] = 0;
    ab_age_male[i]   = 0;
    c_age_female[i]  = 0;
    c_age_male[i]    = 0;
    def_age_female[i]= 0;
    def_age_male[i]  = 0;
  }
  ab_pop_female  = 0;
  ab_pop_male    = 0;
  c_pop_female   = 0;
  c_pop_male     = 0;
  def_pop_female = 0;
  def_pop_male   = 0;
}

// function to calculate proportions of different liver scores in population
void PropLiverScores() {
  
  //ensure accumulators are reset before calculations
  ab_pop_female = 0;
  ab_pop_male   = 0;
  c_pop_female  = 0;
  c_pop_male    = 0;
  def_pop_female = 0;
  def_pop_male   = 0;
  
  for (int i=0; i<na; i++) {
    // Proportion of AB scores (normal liver)
    ab_age_female[i] = female_states(i,0);  //age-specific
    ab_pop_female   += female_states(i,0) * demographic_weights[i];  //overall population
    ab_age_male[i]   =   male_states(i,0);  //age-specific
    ab_pop_male     +=   male_states(i,0) * demographic_weights[i];  //overall population
    // Proportion of C scores (mild fibrosis)
    c_age_female[i] = female_states(i,1);   //age-specific
    c_pop_female   += female_states(i,1) * demographic_weights[i];   //overall population
    c_age_male[i]   =   male_states(i,1);   //age-specific
    c_pop_male     +=   male_states(i,1) * demographic_weights[i];   //overall population
    // Proportion of DEF scores (severe fibrosis)
    def_age_female[i] = female_states(i,2); //age-specific
    def_pop_female   += female_states(i,2) * demographic_weights[i]; //overall population
    def_age_male[i]   =   male_states(i,2); //age-specific
    def_pop_male     +=   male_states(i,2) * demographic_weights[i]; //overall population
  }
  
}


// function to calculate derivatives
List Derivs() {
  
  // Differential equations
  for (int i=0; i<na; i++) {

    if (i == 0) { //first age group
      // AB scores
      dymat_female(i,0) =   eta1 * female_tmp_states(i,1) - tau1 * cur_wb_female[i] * female_tmp_states(i,0) - female_tmp_states(i,0) * da + (female_tmp_states(na-1,0)+female_tmp_states(na-1,1)+female_tmp_states(na-1,2))*da;  //ageing in to balance ageing out from max age of all liver groups
      dymat_male(i,0)   =   eta1 *   male_tmp_states(i,1) - tau1 *   cur_wb_male[i] *   male_tmp_states(i,0) -   male_tmp_states(i,0) * da + (  male_tmp_states(na-1,0)+  male_tmp_states(na-1,1)+  male_tmp_states(na-1,2))*da;  //ageing in to balance ageing out from max age of all liver groups
      // C scores
      dymat_female(i,1) = - eta1 * female_tmp_states(i,1) + tau1 * cur_wb_female[i] * female_tmp_states(i,0) + eta2 * female_tmp_states(i,2) - tau2 * cum_wb_female[i] * female_tmp_states(i,1) - female_tmp_states(i,1) * da; //no ageing into first age category
      dymat_male(i,1)   = - eta1 *   male_tmp_states(i,1) + tau1 *   cur_wb_male[i] *   male_tmp_states(i,0) + eta2 *   male_tmp_states(i,2) - tau2 *   cum_wb_male[i] *   male_tmp_states(i,1) -   male_tmp_states(i,1) * da; //no ageing into first age category
      // DEF scores
      dymat_female(i,2) = - eta2 * female_tmp_states(i,2) + tau2 * cum_wb_female[i] * female_tmp_states(i,1) - female_tmp_states(i,2) * da; //no ageing into first age category
      dymat_male(i,2)   = - eta2 *   male_tmp_states(i,2) + tau2 *   cum_wb_male[i] *   male_tmp_states(i,1) -   male_tmp_states(i,2) * da; //no ageing into first age category
      
    } else {  //subsequent age groups
      // AB scores
      dymat_female(i,0) =   eta1 * female_tmp_states(i,1) - tau1 * cur_wb_female[i] * female_tmp_states(i,0) - female_tmp_states(i,0) * da + female_tmp_states(i-1,0) * da;
      dymat_male(i,0)   =   eta1 *   male_tmp_states(i,1) - tau1 *   cur_wb_male[i] *   male_tmp_states(i,0) -   male_tmp_states(i,0) * da +   male_tmp_states(i-1,0) * da;
      // C scores
      dymat_female(i,1) = - eta1 * female_tmp_states(i,1) + tau1 * cur_wb_female[i] * female_tmp_states(i,0) + eta2 * female_tmp_states(i,2) - tau2 * cum_wb_female[i] * female_tmp_states(i,1) - female_tmp_states(i,1) * da + female_tmp_states(i-1,1)*da;
      dymat_male(i,1)   = - eta1 *   male_tmp_states(i,1) + tau1 *   cur_wb_male[i] *   male_tmp_states(i,0) + eta2 *   male_tmp_states(i,2) - tau2 *   cum_wb_male[i] *   male_tmp_states(i,1) -   male_tmp_states(i,1) * da +   male_tmp_states(i-1,1)*da;
      // DEF scores
      dymat_female(i,2) = - eta2 * female_tmp_states(i,2) + tau2 * cum_wb_female[i] * female_tmp_states(i,1) - female_tmp_states(i,2) * da + female_tmp_states(i-1,2)*da;
      dymat_male(i,2)   = - eta2 *   male_tmp_states(i,2) + tau2 *   cum_wb_male[i] *   male_tmp_states(i,1) -   male_tmp_states(i,2) * da +   male_tmp_states(i-1,2)*da;
    }

  }
  
  // return derivatives
  List derivatives;
  derivatives["dymat_female"] = dymat_female;
  derivatives["dymat_male"]   = dymat_male;
  return derivatives;
  
}


// function to implement RK4 integration
List RK4(NumericMatrix male_states, NumericMatrix female_states, double stepsize){
  
  // set current states in a temporary matrix (use Rcpp::clone() to ensure changes in tmp states don't affect states)
  female_tmp_states = Rcpp::clone(female_states);
  male_tmp_states   = Rcpp::clone(male_states);
  
  // calculate derivatives for temporary states & store
  List derivs_result1 = Derivs();
  k1_female = Rcpp::clone(as<NumericMatrix>(derivs_result1["dymat_female"]));
  k1_male   = Rcpp::clone(as<NumericMatrix>(derivs_result1["dymat_male"]));
  
  for (int i = 0; i < na; i++) {
    for (int j = 0; j < nlg; j++) {
      // move temporary matrices on by stepsize/2
      female_tmp_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/2)*k1_female(i,j);
      male_tmp_states(i,j)   =   male_states(i,j) + (static_cast<double>(stepsize)/2)*k1_male(i,j);
    }
  }
  // calculate derivatives for temporary states and store
  List derivs_result2 = Derivs();
  k2_female = Rcpp::clone(as<NumericMatrix>(derivs_result2["dymat_female"]));
  k2_male   = Rcpp::clone(as<NumericMatrix>(derivs_result2["dymat_male"]));
  
  for (int i = 0; i < na; i++) {
    for (int j = 0; j < nlg; j++) {
      // move temporary matrices on by stepsize/2
      female_tmp_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/2)*k2_female(i,j);
      male_tmp_states(i,j)   =   male_states(i,j) + (static_cast<double>(stepsize)/2)*k2_male(i,j);
    }
  }
  // calculate derivatives for temporary states & store
  List derivs_result3 = Derivs();
  k3_female = Rcpp::clone(as<NumericMatrix>(derivs_result3["dymat_female"]));
  k3_male   = Rcpp::clone(as<NumericMatrix>(derivs_result3["dymat_male"]));
  
  for (int i = 0; i < na; i++) {
    for (int j = 0; j < nlg; j++) {
      // move temporary matrices on by stepsize
      female_tmp_states(i,j) = female_states(i,j) + stepsize*k3_female(i,j);
      male_tmp_states(i,j)   =   male_states(i,j) + stepsize*k3_male(i,j);
    }
  }
  // calculate derivatives for temporary states & store
  List derivs_result4 = Derivs();
  k4_female = Rcpp::clone(as<NumericMatrix>(derivs_result4["dymat_female"]));
  k4_male   = Rcpp::clone(as<NumericMatrix>(derivs_result4["dymat_male"]));
  
  //update the state matrices
  for (int i = 0; i < na; i++) {
    for (int j = 0; j < nlg; j++) {
      female_updated_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/6)*( k1_female(i,j) + 2*k2_female(i,j) + 2*k3_female(i,j) + k4_female(i,j) );
      male_updated_states(i,j)   =   male_states(i,j) + (static_cast<double>(stepsize)/6)*( k1_male(i,j)   + 2*k2_male(i,j)   + 2*k3_male(i,j)   + k4_male(i,j) );
    }
  }
  
  // Create a List to store both updated states and return it
  List result;
  result["male_updated_states"]   = male_updated_states;
  result["female_updated_states"] = female_updated_states;
  return result;
  
}


// Function to explore relationship of AB <--> C & C <--> DEF
void DerivsRatios() {
  
  // ensure accumulators set to 0
  AB_C_ratio_pop_female  = 0;
  C_DEF_ratio_pop_female = 0;
  
  for (int i = 0; i < na; i++) {
    
    // FEMALES -----------------------------------------------------------------
    // AB --> C
    AB_C_female[i] = tau1 * cur_wb_female[i] * female_tmp_states(i,0);
    // C --> AB
    C_AB_female[i] = eta1 * female_tmp_states(i,1);
    // Ratio of (AB --> C) / (C --> AB)
    AB_C_ratio_female[i] = AB_C_female[i] / C_AB_female[i];
    // Population mean for given timepoint
    AB_C_ratio_pop_female += AB_C_ratio_female[i] * demographic_weights[i];
    
    // C --> DEF
    C_DEF_female[i] =  tau2 * cum_wb_female[i] * female_tmp_states(i,1);
    // DEF --> C
    DEF_C_female[i] =  eta2 * female_tmp_states(i,2);
    // Ratio of (C --> DEF) / (DEF --> C)
    C_DEF_ratio_female[i] = C_DEF_female[i] / DEF_C_female[i];
    // Population mean for given timepoint
    C_DEF_ratio_pop_female += C_DEF_ratio_female[i] * demographic_weights[i];
    
  }
  
}


// [[Rcpp::export]]
List RunMorbidityModel(NumericVector pars, int runtime, double stepsize, NumericVector wormburden_female, NumericVector wormburden_male,
                       NumericVector cumulative_wormburden_female, NumericVector cumulative_wormburden_male) {
  // List RunMorbidityModel(NumericVector pars, int runtime, double stepsize, NumericMatrix wormburden_female, NumericMatrix wormburden_male, 
  //                        NumericMatrix cumulative_wormburden_female, NumericMatrix cumulative_wormburden_male) {
  
  //define parameters to be varied
  tau1 = exp(pars[0]);
  tau2 = exp(pars[1]);
  eta1 = exp(pars[2]);
  eta2 = exp(pars[3]);
  
  // calculate no. time steps
  nsteps = static_cast<int>(runtime/stepsize);  //force nsteps to be integer
  
  // define storage vectors for output (after no. time steps is defined)
  NumericVector time_out(nsteps);
  NumericVector ab_pop_female_out(nsteps), ab_pop_male_out(nsteps), c_pop_female_out(nsteps), c_pop_male_out(nsteps), def_pop_female_out(nsteps), def_pop_male_out(nsteps);
  NumericMatrix ab_age_female_out(nsteps, na), ab_age_male_out(nsteps, na), c_age_female_out(nsteps, na), c_age_male_out(nsteps, na), def_age_female_out(nsteps, na), def_age_male_out(nsteps, na);

  NumericVector AB_C_ratio_pop_female_out(nsteps), C_DEF_ratio_pop_female_out(nsteps);
  
  // run initialization functions
  AgeMidpoints();          //age midpoints
  DemographicWeights();    //demographic weights
  InitialiseStates();      //state variables
  InitialiseLiverScores(); //liver score variables
  
  // store initial matrix for output
  for (int i = 0; i < na; i++) {    //loop over age groups
    for (int j = 0; j < nlg; j++) { //loop over liver groups
      female_init_states(i,j) = female_states(i,j);
      male_init_states(i,j)   =   male_states(i,j);
    }
  }
  
  // calculate mean proportion of different liver scores based on initial states
  PropLiverScores();
  
  // store output for time 0
  time_out[0]      = 0;
  ab_pop_female_out[0]    = ab_pop_female;
  ab_pop_male_out[0]      = ab_pop_male;
  c_pop_female_out[0]     = c_pop_female;
  c_pop_male_out[0]       = c_pop_male;
  def_pop_female_out[0]   = def_pop_female;
  def_pop_male_out[0]     = def_pop_male;
  ab_age_female_out(0,_)  = ab_age_female;
  ab_age_male_out(0,_)    = ab_age_male;
  c_age_female_out(0,_)   = c_age_female;
  c_age_male_out(0,_)     = c_age_male;
  def_age_female_out(0,_) = def_age_female;
  def_age_male_out(0,_)   = def_age_male;
  
  // set time-fixed worm burdens
  cur_wb_female = Rcpp::clone(wormburden_female);
  cur_wb_male   = Rcpp::clone(wormburden_male);
  cum_wb_female = Rcpp::clone(cumulative_wormburden_female);
  cum_wb_male   = Rcpp::clone(cumulative_wormburden_male);
  
  // loop through time steps
  for (int h = 1; h < nsteps; h++){
    
    // update time
    time_out[h] = h*stepsize;
    
    // //extract time-specific current & cumulative worm burdens from input from transmission model
    // cur_wb_female = wormburden_female(h,_);
    // cur_wb_male   = wormburden_male(h,_);
    // cum_wb_female = cumulative_wormburden_female(h,_);
    // cum_wb_male   = cumulative_wormburden_male(h,_);

    // perform RK4 integration
    List updated_states   = RK4(male_states, female_states, stepsize);
    male_updated_states   = as<NumericMatrix>(updated_states["male_updated_states"]);
    female_updated_states = as<NumericMatrix>(updated_states["female_updated_states"]);
    
    // replace state matrices with updated matrices (done in long form)
    for (int i = 0; i < na; i++) {
      for (int j = 0; j < nlg; j++) {
        female_states(i,j) = female_updated_states(i,j);
        male_states(i,j)   =   male_updated_states(i,j);
      }
    }

    // calculate updated proportions of liver scores
    PropLiverScores();
    
    ab_pop_female_out[h]    = ab_pop_female;
    ab_pop_male_out[h]      = ab_pop_male;
    c_pop_female_out[h]     = c_pop_female;
    c_pop_male_out[h]       = c_pop_male;
    def_pop_female_out[h]   = def_pop_female;
    def_pop_male_out[h]     = def_pop_male;
    ab_age_female_out(h,_)  = ab_age_female;
    ab_age_male_out(h,_)    = ab_age_male;
    c_age_female_out(h,_)   = c_age_female;
    c_age_male_out(h,_)     = c_age_male;
    def_age_female_out(h,_) = def_age_female;
    def_age_male_out(h,_)   = def_age_male;
    
    // calculate ratios 
    // DerivsRatios();
    // AB_C_ratio_pop_female_out[h]  = AB_C_ratio_pop_female;
    // C_DEF_ratio_pop_female_out[h] = C_DEF_ratio_pop_female;
    
    
  } //end of time loop
  
  
  // return outputs
  return Rcpp::List::create(
    Rcpp::Named("time") = time_out,
    Rcpp::Named("age") = age,
    Rcpp::Named("female_init_states") = female_init_states, 
    Rcpp::Named("male_init_states")   = male_init_states, 
    Rcpp::Named("ab_pop_female")      = ab_pop_female_out, 
    Rcpp::Named("ab_pop_male")        = ab_pop_male_out, 
    Rcpp::Named("c_pop_female")       = c_pop_female_out, 
    Rcpp::Named("c_pop_male")         = c_pop_male_out, 
    Rcpp::Named("def_pop_female")     = def_pop_female_out, 
    Rcpp::Named("def_pop_male")       = def_pop_male_out, 
    Rcpp::Named("ab_age_female")      = ab_age_female_out, 
    Rcpp::Named("ab_age_male")        = ab_age_male_out, 
    Rcpp::Named("c_age_female")       = c_age_female_out, 
    Rcpp::Named("c_age_male")         = c_age_male_out, 
    Rcpp::Named("def_age_female")     = def_age_female_out,
    Rcpp::Named("def_age_male")       = def_age_male_out//,
    // Rcpp::Named("AB_C_ratio_female")  = AB_C_ratio_pop_female_out,
    // Rcpp::Named("C_DEF_ratio_female") = C_DEF_ratio_pop_female_out
    
    );
  
}
