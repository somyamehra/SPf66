#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

/* Parameters:
 *    N_OBS: length of study period (in windows)
 *    N_AGE: age at onset of study period (in windows)
 *    TIME_STEP: length of each window
 *    ETA: hypnozoite clearance rate
 *    NU: average number of relapses per bite
 *    PCLIN: prob that each hypnozoite activation/reinfection event during the 
             study period gives rise to clinical symptoms
 *    PREL: prob that a sporozoite forms a hypnozoite (c.f immediate development)
 *    FORI: vector of length N_OBS+N_WINDOWS for time-dependent FORI 
 *    THRESHOLD_AGE: model age-dependence on FORI as a step function  
 *    U: proportion reduction in FORI before THRESHOLD_AGE
*/
 

/* Evaluate the log of the joint PGF [9] of infection states over the study
 * period for all row-vectors of a given matrix
 * 
 * Additional arguments:
 *    Z: matrix with N_OBS columns
*/

// [[Rcpp::export]]
vec log_joint_PGF 
  (mat Z, int N_OBS, int N_AGE, double TIME_STEP, double ETA, double NU, 
   double PCLIN, double PREL, vec FORI) {


  if (N_OBS <= 0 | N_AGE <=0 | TIME_STEP <=0 | ETA <= 0 | NU <= 0 | PCLIN <= 0 | PREL <=0) {
    Rcpp::stop("parameters must be non-negative");
  } 
  
  if (PCLIN > 1 | PREL > 1) {
    Rcpp::stop("PCLIN, PREL cannot be greater than 1");
  } 
  
  
  if (FORI.n_elem != N_AGE + N_OBS) {
    Rcpp::stop("FORI vector must have length N_AGE+N_OBS");
  }
  
  if (Z.n_cols != N_OBS) {
    Rcpp::stop("matrix Z must have N_OBS columns");
  }
  
  vec exp_eta_study = exp(-ETA*TIME_STEP*linspace(1, N_OBS, N_OBS));
  vec exp_eta_age_diff = exp(-ETA*TIME_STEP*linspace(N_AGE-1, 0, N_AGE));
  
  // contribution from bites before study period
  vec z_sum_pre_study =  exp(-ETA*TIME_STEP) - exp(-ETA*TIME_STEP*(N_OBS+1)) - 
    (1-exp(-ETA*TIME_STEP))*sum(Z.each_row() % exp_eta_study.t(), 1);
  
  mat p1 = log1p(-(1-exp(-ETA*TIME_STEP))/(1 + NU*PCLIN*PREL*(z_sum_pre_study * exp_eta_age_diff.t())));
  
  p1.each_row() %= FORI.subvec(0, N_AGE-1).t() / ETA;
  
  // contribution from bites in study period
  
  mat z_sum_study = -NU*PREL*PCLIN*(1-exp(-ETA*TIME_STEP))*reverse(cumsum(reverse(Z.each_row() % exp_eta_study.t(), 1), 1), 1);
  
  z_sum_study.each_row() %= pow(exp_eta_study, -1).t();
  
  z_sum_study.each_row() += 1 + NU*PREL*PCLIN*(1-reverse(exp_eta_study.t()));
  
  mat z_function_1 = (1 - PCLIN + PCLIN*Z)/(ETA*(1 + NU*PCLIN*PREL - NU*PCLIN*PREL*Z));
  mat coeff_study_1 = z_function_1.each_row() % FORI.subvec(N_AGE, N_OBS+N_AGE-1).t();

  mat z_function_2 = (PCLIN - PCLIN*Z)/(ETA*(1 + NU*(1-PREL) + NU*PCLIN*PREL - NU*PCLIN*PREL*Z));
  mat coeff_study_2 = z_function_2.each_row() % FORI.subvec(N_AGE, N_OBS+N_AGE-1).t();
  
  mat p2= log1p(-(1+NU*PCLIN*PREL-NU*PCLIN*PREL*Z)*(1-exp(-ETA*TIME_STEP))/z_sum_study) % coeff_study_1;
  mat p3= log1p(-(1+NU*(1-PREL)+NU*PCLIN*PREL-NU*PCLIN*PREL*Z)*(1-exp(-ETA*TIME_STEP))/(NU*(1-PREL)+z_sum_study)) % coeff_study_2;
  
  return(-(sum(p1, 1) + sum(p2, 1) + sum(p3, 1) + TIME_STEP*sum(FORI)));
  
}

/* Calculate the likelihood of a binary sequence c by applying the 
 * inclusion-exclusion principle to the joint PGF [9].
 * Rather than using c directly, we pre-evaluate the sets S_+(c) (Eq [7]) 
 * and S_-(c) (Eq [8]) and use Eq [6].
 * 
 * Additional arguments: for a given binary sequence c
 *    S_plus: row-vectors correspond to each element of S_+(c) 
 *    S_minus: row-vectors correspond to each element of S_-(c) 
 */

// [[Rcpp::export]]
double inc_exc_likelihood 
  (mat S_plus, mat S_minus, int N_AGE, double PCLIN, vec FORI_full,
   double ETA, double NU, double PREL, int N_OBS, double TIME_STEP,
   int THRESHOLD_AGE, double U) {
  
  vec FORI = FORI_full.tail(N_AGE+N_OBS);
  
  FORI.subvec(0, THRESHOLD_AGE-1) *= U;
  
  if (S_minus.n_rows==0) {
    return(exp(log_joint_PGF(S_plus, N_OBS, N_AGE, TIME_STEP, ETA, NU, PCLIN, PREL, FORI).front()));
  }
  
  vec f_plus = log_joint_PGF(S_plus, N_OBS, N_AGE, TIME_STEP, ETA, NU, PCLIN, PREL, FORI);
  
  vec f_minus = log_joint_PGF(S_minus, N_OBS, N_AGE, TIME_STEP, ETA, NU, PCLIN, PREL, FORI);
  
  double max_val = std::max(f_plus.max(), f_minus.max());
  
  f_plus = exp(f_plus - max_val);
  f_minus = exp(f_minus - max_val);
  
  // |S_+(c)| - |S_-(c)| = |c| mod 2
  if (f_minus.n_elem == f_plus.n_elem) {
    return(exp(max_val) * sum(f_plus - f_minus));
  } else {
    return(exp(max_val)*(sum(f_plus.subvec(1, f_minus.n_elem) - f_minus) + f_plus.front()));
  }
}

