library(parallel)
library(Rcpp)
library(RcppArmadillo)

sourceCpp('Model_calibration_code/clin_inf_likelihood.cpp')

# ============================================================================
# Generate a set of S vectors given a list of ternary infection states

generate_S_vectors <- function(inf_states_by_VIN) {
  S_vectors <- lapply(inf_states_by_VIN, function(inf_state) {
    obs_group <- rep(1, length(inf_state))
    curr_group <- 1
    for (i in 2:length(inf_state)) {
      curr_group <- curr_group + as.numeric(inf_state[i-1]!="B" | inf_state[i]!="B" )
      obs_group[i] <- curr_group
    }
    bunched_inf_state <- inf_state[!duplicated(obs_group)]
    x <- rep(list(0), curr_group)
    x[bunched_inf_state=="M"] <- list(NA)
    x[bunched_inf_state %in% c("B", "C")] <- list(c(0, 1))
    z_vals <- expand.grid(x)
    pos_vals <- (rowSums(z_vals, na.rm = T)%%2-sum(bunched_inf_state %in% c("B", "C"))%%2==0)
    z_vals[is.na(z_vals)] <- 1
    return(list(plus=as.matrix(z_vals[pos_vals, obs_group]),
                minus=as.matrix(z_vals[!pos_vals, obs_group])))
  })
  return(list(S_plus=lapply(S_vectors, function(x) as.matrix(x[["plus"]])),
              S_minus=lapply(S_vectors, function(x) as.matrix(x[["minus"]]))))
}

# ============================================================================
# Calculate likelihood for a given cohort

cohort_likelihood <- function(S_plus, S_minus, N_AGES, PCLIN, FORI, ETA, NU, PREL,
                              N_OBS, TIME_STEP, THRESHOLD_AGE, U, N_CORES_PER_CHAIN) {
  
  indiv_likelihoods <- mcmapply(inc_exc_likelihood,
                                S_plus,
                                S_minus,
                                N_AGES,
                                PCLIN,
                                MoreArgs = list(FORI=FORI, ETA=ETA, NU=NU, PREL=PREL, 
                                                N_OBS=N_OBS, TIME_STEP=TIME_STEP,
                                                THRESHOLD_AGE=THRESHOLD_AGE, U=U),
                                mc.cores = N_CORES_PER_CHAIN)
  
  # numerical instability for implausible params: adjust negative/zero likelihoods
  indiv_likelihoods[indiv_likelihoods<=0] <- 1e-16
  
  return(sum(log10(indiv_likelihoods)))
}

# ============================================================================

# Calculate antidisease masking probability
calculate_pclin <- function(N_AGES, LOGIT_RHO, LOG_GAMMA, MIN_AGE, MAX_AGE) {
  1 - (1-LaplacesDemon::invlogit(LOGIT_RHO))*
    ((N_AGES-MIN_AGE)/(MAX_AGE-MIN_AGE))^exp(LOG_GAMMA)
}

# Calculate time-dependent FORI
time_dependent_fori <- function(LAMBDA1, LAMBDA2, SHIFT_WINDOW, SEASONALITY) {
  c(head(SEASONALITY, -(N_OBS-SHIFT_WINDOW))*LAMBDA1,
    tail(SEASONALITY, N_OBS-SHIFT_WINDOW)*LAMBDA2)
}