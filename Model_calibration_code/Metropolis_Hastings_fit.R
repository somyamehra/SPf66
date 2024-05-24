library(parallel)
library(dplyr)
library(tidyr)

# ==============================================================================
# Metropolis-Hastings algorithm for model fitting
# nested parallelistion: parallelisation over chains, calculation of
# cohort likelihood within each chain
Metropolis_Hastings <- function(inf_states_by_VIN, N_AGES, PREL, N_OBS, TIME_STEP,
                                SEASONALITY, SHIFT_WINDOW, THRESHOLD_AGE, U,
                                N_CHAIN, N_ITER, N_CORES_PER_CHAIN,
                                LAMBDA_PROP_SD, NU_PROP_SD, ETA_PROP_SD,
                                LOGIT_RHO_PROP_SD, LOG_GAMMA_PROP_SD,
                                LOG_GAMMA_PRIOR_SD, LOGIT_RHO_PRIOR_SD) {
  
  source("Model_calibration_code/clin_inf_likelihood_helper_func.R")
  
  # likelihood computed using the inclusion-exclusion principle
  # generate relevant sets of S-vectors for each individual  
  S_vectors <- generate_S_vectors(inf_states_by_VIN)
  S_plus <- S_vectors[["S_plus"]]
  S_minus <- S_vectors[["S_minus"]]
  
  # https://stackoverflow.com/questions/8413188/can-i-nest-parallelparlapply
  # parallelise over chains
  # this works on Mac but a different workaround might be needed for Windoes 
  cl <- makeCluster(N_CHAINS)
  
  # export objects from local function environment
  clusterExport(cl, c('S_plus', 'S_minus', 'N_AGES', 'PREL', 'N_OBS', 'TIME_STEP',
                      'SEASONALITY', 'SHIFT_WINDOW', 'THRESHOLD_AGE', 'U',
                      'N_ITER', 'N_CORES_PER_CHAIN',
                      'LAMBDA_PROP_SD', 'NU_PROP_SD', 'ETA_PROP_SD',
                      'LOGIT_RHO_PROP_SD', 'LOG_GAMMA_PROP_SD',
                      'LOG_GAMMA_PRIOR_SD', 'LOGIT_RHO_PRIOR_SD'),
                envir=environment())
  
  MCMC_results <- parLapply(cl, 1:length(cl), function(CHAIN) {
    
    source("Model_calibration_code/clin_inf_likelihood_helper_func.R")
    
    MIN_AGE <- min(N_AGES)
    MAX_AGE <- max(N_AGES)
    
    param_vals <- data.frame(matrix(ncol = 9, nrow = N_ITER))
    colnames(param_vals) <- c("CHAIN", "ITER", "LAMBDA1", "LAMBDA2", "NU", "ETA", "LOGIT_RHO", "LOG_GAMMA", "LOG_LIK")
    param_vals[, "CHAIN"] <- CHAIN
    param_vals[, "ITER"] <- 1:N_ITER
    
    # initialise parameter values
    param_vals[1, "LAMBDA1"] <-  runif(1, 0.1/365, 1/365)
    param_vals[1, "LAMBDA2"] <-  runif(1, 0.1/365, 1/365)
    param_vals[1, "NU"] <-  runif(1, 0.5, 8)
    param_vals[1, "ETA"] <- runif(1, 1/500, 1/50)
    param_vals[1, "LOGIT_RHO"] <- rnorm(1, sd=LOGIT_RHO_PROP_SD)
    param_vals[1, "LOG_GAMMA"] <- rnorm(1, sd=LOG_GAMMA_PROP_SD)
    param_vals[1, "LOG_LIK"] <- 
      cohort_likelihood(S_plus, S_minus, N_AGES, 
                        calculate_pclin(N_AGES, param_vals[1, "LOGIT_RHO"], 
                                        param_vals[1, "LOG_GAMMA"], MIN_AGE, MAX_AGE),
                        time_dependent_fori(param_vals[1, "LAMBDA1"], param_vals[1, "LAMBDA2"],
                                            SHIFT_WINDOW, SEASONALITY),
                        param_vals[1, "ETA"], param_vals[1, "NU"], PREL, 
                        N_OBS, TIME_STEP, THRESHOLD_AGE, U, N_CORES_PER_CHAIN)
    
    for (i in 2:N_ITER) {
      # sample parameter values from proposal distrbution
      lambda1_prop <- rnorm(1, param_vals[i-1, "LAMBDA1"], LAMBDA_PROP_SD)
      if (lambda1_prop<0) lambda1_prop <- param_vals[i-1, "LAMBDA1"]
      
      lambda2_prop <- rnorm(1, param_vals[i-1, "LAMBDA2"], LAMBDA_PROP_SD)
      if (lambda2_prop<0) lambda2_prop <- param_vals[i-1, "LAMBDA2"]
      
      nu_prop <- rnorm(1, param_vals[i-1, "NU"], NU_PROP_SD)
      if (nu_prop<0) nu_prop <- param_vals[i-1, "NU"]
      
      eta_prop <- rnorm(1, param_vals[i-1, "ETA"], ETA_PROP_SD)
      if (eta_prop<0) eta_prop <- param_vals[i-1, "ETA"]
      
      logit_rho_prop <- rnorm(1, param_vals[i-1, "LOGIT_RHO"], LOGIT_RHO_PROP_SD)
      
      log_gamma_prop <- rnorm(1, param_vals[i-1, "LOG_GAMMA"], LOG_GAMMA_PROP_SD)
      
      # calculate likelihood ratio
      likelihood_prop <- 
        cohort_likelihood(S_plus, S_minus, N_AGES, 
                          calculate_pclin(N_AGES, logit_rho_prop, log_gamma_prop, 
                                          MIN_AGE, MAX_AGE),
                          time_dependent_fori(lambda1_prop, lambda2_prop,
                                              SHIFT_WINDOW, SEASONALITY),
                          eta_prop, nu_prop, PREL, N_OBS, TIME_STEP, 
                          THRESHOLD_AGE, U, N_CORES_PER_CHAIN)
      
      likelihood_ratio <- exp(likelihood_prop-param_vals[i-1, "LOG_LIK"])*
        dnorm(log_gamma_prop, sd = LOG_GAMMA_PRIOR_SD)*
        dnorm(logit_rho_prop, sd = LOGIT_RHO_PRIOR_SD)/
        (dnorm(param_vals[i-1, "LOG_GAMMA"], sd = LOG_GAMMA_PRIOR_SD)*
           dnorm(param_vals[i-1, "LOGIT_RHO"], sd = LOGIT_RHO_PRIOR_SD))
      
      # accept or reject
      if (runif(1)<=likelihood_ratio) {
        param_vals[i, "LAMBDA1"] <- lambda1_prop
        param_vals[i, "LAMBDA2"] <- lambda2_prop
        param_vals[i, "NU"] <- nu_prop
        param_vals[i, "ETA"] <- eta_prop
        param_vals[i, "LOG_LIK"] <- likelihood_prop
        param_vals[i, "LOGIT_RHO"] <- logit_rho_prop
        param_vals[i, "LOG_GAMMA"] <- log_gamma_prop
      } else {
        param_vals[i, 3:9] <- param_vals[i-1, 3:9]
      }
    }
    
    return(param_vals)
    
  })
  
  return(dplyr::bind_rows(MCMC_results))
}

# ==============================================================================
# calculate Gelman-Rubin convergence diagnostic

calculate_Gelman_Rubin <- function(MCMC_results, N_ITER, BURNIN, N_CHAINS, PARAMS) {
  
  MCMC_results <- MCMC_results %>% subset(ITER>=BURNIN)
  
  within_chain_var <- MCMC_results %>% dplyr::select(c("CHAIN", PARAMS)) %>%
    group_by(CHAIN) %>% summarise_all(var) %>% dplyr::select(-CHAIN)
  
  mean_within_chain_var <- colMeans(within_chain_var)
  
  within_chain_mean <- MCMC_results %>% dplyr::select(c("CHAIN", PARAMS)) %>%
    group_by(CHAIN) %>% summarise_all(mean) %>% dplyr::select(-CHAIN)
  
  var_within_chain_mean <- matrixStats::colVars(as.matrix(within_chain_mean))
  
  return(sqrt(1 - 1/N_ITER + (1+1/N_CHAINS)*var_within_chain_mean/mean_within_chain_var))
  
}
