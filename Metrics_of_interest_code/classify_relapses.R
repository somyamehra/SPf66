library(Rcpp)
library(RcppArmadillo)

source('Model_calibration_code/clin_inf_likelihood_helper_func.R')
sourceCpp('Metrics_of_interest_code/clin_inf_primary_likelihood.cpp')


# probabilistically classify each symptomatic episode as a relapse (defined to be the 
# probability that it is caused by hypnozoite activation event(s) only)

classify_relapses <- function(inf_states_by_VIN, N_AGES, PCLIN, FORI, ETA, NU, PREL) {
  S_vectors <- generate_S_vectors(inf_states_by_VIN)
  S_plus <- S_vectors[["S_plus"]]
  S_minus <- S_vectors[["S_minus"]]
  
  
  keep_VIN <- names(Filter(function(x) sum(x %in% c("C", "B"))>0, inf_states_by_VIN))
  
  # identify windows to perform classification 
  # may be aggregated due to prophylactic bunching
  TARGET_WINDOWS <- lapply(inf_states_by_VIN[keep_VIN], function(inf_state) {
    obs_group <- rep(1, length(inf_state))
    curr_group <- 1
    for (i in 2:length(inf_state)) {
      curr_group <- curr_group + as.numeric(inf_state[i-1]!="B" | inf_state[i]!="B" )
      obs_group[i] <- curr_group
    }
    inf_groups <- split(which(inf_state %in% c("C", "B")), 
                        obs_group[inf_state %in% c("C", "B")])
    names(inf_groups) <- sapply(inf_groups, max)
    return(inf_groups)
  })
  
  # probabilistic classification
  prob_relapse_only <- lapply(keep_VIN, function(my_VIN) {
    sapply(TARGET_WINDOWS[[my_VIN]], function(WINDOWS) {
      no_prim_prob <- inc_exc_likelihood_noprim_window(
        S_plus[[my_VIN]], S_minus[[my_VIN]], N_AGES[my_VIN], PCLIN[my_VIN],
        PREL=PREL, N_OBS=N_OBS, TIME_STEP=TIME_STEP, ETA=ETA, NU=NU, 
        FORI=FORI, NO_PRIM_WINDOW=WINDOWS)
      
      recurrence_prob <- inc_exc_likelihood(
        S_plus[[my_VIN]], S_minus[[my_VIN]], N_AGES[my_VIN], PCLIN[my_VIN],
        PREL=PREL, N_OBS=N_OBS, TIME_STEP=TIME_STEP, ETA=ETA, NU=NU, 
        FORI=FORI, THRESHOLD_AGE=1, U=1)
      return(no_prim_prob/recurrence_prob)})
  }) 
  names(prob_relapse_only) <- keep_VIN
  
  recurrence_classification <- prob_relapse_only %>% 
    lapply(function(x) data.frame(WINDOW=names(x), PROB_RELAPSE=x)) %>% 
    bind_rows(.id="VIN") %>% 
    mutate(WINDOW=as.numeric(WINDOW),
           PROB_RELAPSE=ifelse(PROB_RELAPSE<=1 & PROB_RELAPSE>=0, PROB_RELAPSE, NA)) 
  
  
  return(recurrence_classification)
} 
