# Metrics of epidemiological interest
# Assumption: hypnozoite reservoir has reached stationarity under a constant FOI

# ==============================================================================
# Size of the hypnozoite reservoir

hyp_reservoir_stationary_pmf <- function(HMAX, ETA, NU, PREL, LAMBDA) {
  dnbinom(0:HMAX, mu=NU*PREL*LAMBDA/ETA, size = LAMBDA/ETA)
}

hyp_reservoir_stationary_cdf <- function(HMAX, ETA, NU, PREL, LAMBDA) {
  pnbinom(0:HMAX, mu=NU*PREL*LAMBDA/ETA, size = LAMBDA/ETA)
}

# ==============================================================================
# Time to spontaneous clearance of hypnozoite reservoir
# Assuming mosquito-to-human transmission temporarily interrupted

hyp_reservoir_stationary_clearance_cdf <- function(t, ETA, NU, PREL, LAMBDA) {
  (1+NU*PREL*exp(-ETA*t))^-(LAMBDA/ETA)
}

# years to clear hypnozoite reservoir with probability THRESHOLD
time_to_clear_hyp <- function(THRESHOLD, ETA, NU, PREL, LAMBDA) {
  max(0, -log((THRESHOLD^(-ETA/LAMBDA)-1)/(NU*PREL))/(ETA*365))
}

# ==============================================================================
# Accuracy of recent recurrence (within t days) as a predictor of hyp carriage

recent_recur_accuracy_stationary <- function(t, ETA, NU, PREL, LAMBDA) {
  c(Specificity=exp(-LAMBDA*t*NU/(1+NU)),
    Hyp_Carriage=1-(1+NU*PREL)^(-LAMBDA/ETA),
    Recent_Recurrence=1-exp(-LAMBDA*t)*(1+NU*PREL*(1-exp(-ETA*t)))^(-LAMBDA/ETA)*
      (1-(1+NU)*(1-exp(-ETA*t))/(1+NU-NU*PREL*exp(-ETA*t)))^(-LAMBDA/(ETA*(1+NU))),
    Sensitivity=(1-((1+NU*PREL)^(-LAMBDA/ETA))*(1-exp(-LAMBDA*t*NU/(1+NU)))-
                   exp(-LAMBDA*t)*(1+NU*PREL*(1-exp(-ETA*t)))^(-LAMBDA/ETA)*
                   (1-(1+NU)*(1-exp(-ETA*t))/
                      (1+NU-NU*PREL*exp(-ETA*t)))^(-LAMBDA/(ETA*(1+NU))))/
      (1-(1+NU*PREL)^(-LAMBDA/ETA))) }

# ==============================================================================
# Specificity of recent recurrence (within t days) as a predictor of hyp carriage,
# given the force of inoculation is sampled randomly from a Gamma distribution 
# with mean LAMBDA_MEAN and scale parameter THETA

recent_recur_specificity_stationary_het <- function(t, ETA, NU, PREL, LAMBDA_MEAN, THETA) {
  (1+THETA*t*NU/(1+NU))^(-LAMBDA_MEAN/THETA)
}

# Proportion of bites experienced by the 20% of individuals subject to the highest
# transmission intensity given the force of inoculation in the population is 
# Gamma-distributed with population mean LAMBDA_MEAN and scale parameter THETA
prop_bites_top_20 <- function(THETA, LAMBDA_MEAN) {
  0.2 + dgamma(x = qgamma(0.8, shape=LAMBDA_MEAN/THETA, scale=THETA),
               shape = LAMBDA_MEAN/THETA+1, scale=THETA)*THETA
}

  
# ==============================================================================
# Probability that first recurrence occurs within observation window n  
time_to_first_recur <- function(n, TIME_STEP, LAMBDA, NU, PREL, ETA) {
  exp(-LAMBDA*n*TIME_STEP)*(1+NU*PREL*(1-exp(-ETA*n*TIME_STEP)))^(-LAMBDA/ETA)*
    (1-(1+NU)*(1-exp(-ETA*n*TIME_STEP))/
       (1+NU-NU*PREL*exp(-ETA*n*TIME_STEP)))^(-LAMBDA/(ETA*(1+NU)))
}

# Probability of a baseline recurrence in window 1, subsequent recurrence in window n
next_recurrence_after_baseline <- function(n, TIME_STEP, LAMBDA, NU, PREL, ETA) {
  time_to_first_recur(n-2, TIME_STEP, LAMBDA, NU, PREL, ETA)-
    2*time_to_first_recur(n-1, TIME_STEP, LAMBDA, NU, PREL, ETA) +
    time_to_first_recur(n, TIME_STEP, LAMBDA, NU, PREL, ETA)
}

# Conditional on a baseline recurrence in window 1, probability that the subsequent
# recurrence occurs in window n AND is derived from the same batch
prob_same_batch <- function(n, TIME_STEP, LAMBDA, NU, PREL, ETA) {
  (-1/ETA*(log(1+NU*PREL*(1-exp(-ETA*(n-2)*TIME_STEP))) + 
             log(1+NU*PREL*(1-exp(-ETA*n*TIME_STEP))) -
             2*log(1+NU*PREL*(1-exp(-ETA*(n-1)*TIME_STEP)))) +
     -1/(ETA*(1+NU))*(log(1+NU-NU*PREL*exp(-ETA*(n-2)*TIME_STEP)) + 
                        log(1+NU-NU*PREL*exp(-ETA*n*TIME_STEP))-
                        2*log(1+NU-NU*PREL*exp(-ETA*(n-1)*TIME_STEP))))*
    time_to_first_recur(n-2, TIME_STEP, LAMBDA, NU, PREL, ETA)*LAMBDA
}

# Conditional on a baseline recurrence in window 1 and the subsequent recurrence
# in window n, probability that both recurrences are derived from the same batch
conditional_prob_same_batch <- function(n, TIME_STEP, LAMBDA, NU, PREL, ETA) {
  prob_same_batch(n, TIME_STEP, LAMBDA, NU, PREL, ETA)/
    next_recurrence_after_baseline(n, TIME_STEP, LAMBDA, NU, PREL, ETA)
}

# Conditional on a baseline recurrence in window 1, probability that the next
# recurrence occurs in window n
conditional_prob_next_recurrence <- function(n, TIME_STEP, LAMBDA, NU, PREL, ETA) {
  next_recurrence_after_baseline(n, TIME_STEP, LAMBDA, NU, PREL, ETA)/
    (1-time_to_first_recur(1, TIME_STEP, LAMBDA, NU, PREL, ETA))
}

# ==============================================================================
# FOI required to yield recurrence with probability p_recur in a window of W days
fori_required_for_recur <- function(W, p_recur, NU, PREL, ETA, PCLIN) {
  -log(1-p_recur)/
    (W*(1 - (1-PCLIN)/(1+NU*PREL*PCLIN) - PCLIN/(1+NU*(1-PREL*(1-PCLIN)))) + 
       log(1 + NU*PREL*PCLIN*(1-exp(-ETA*W)))*(1 - (1-PCLIN)/(1+NU*PREL*PCLIN))/ETA +
       log((1+NU*(1-PREL))/(1+NU*(1-PREL)+NU*PREL*PCLIN*(1-exp(-ETA*W))))*
       PCLIN/(ETA*(1+NU*(1-PREL*(1-PCLIN)))))
}

