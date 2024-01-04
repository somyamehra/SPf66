library(dplyr)

# simulate data under within-host model

# generate bite times according to a non-homogeneous Poisson process
simulate_bite_times <- function(LAMBDA, TIME_STEP, MAX_TIME) {
  MAX_LAMBDA <- max(LAMBDA) 
  max_bites <- rpois(1, MAX_LAMBDA*TIME_STEP*MAX_TIME)
  bite_times <- runif(max_bites, max=MAX_TIME)
  
  if (max_bites>0) {
    keep_thinned <- runif(max_bites) < LAMBDA[ceiling(bite_times)]/MAX_LAMBDA
    bite_times <- bite_times[keep_thinned]*TIME_STEP
  }
  
  return(bite_times)
}

# simulate bite times associated with each relapse
# R is the size parameter for the negative binomial distribution
simulate_recurrence_times <- function(NU, R, ETA, PREL) {
  n_sporozoites <- rnbinom(1, mu=NU, size=R)
  n_hypnozoites <- rbinom(1, n_sporozoites, PREL)
  
  recur_times <- rexp(n_hypnozoites, ETA)
  if (n_hypnozoites<n_sporozoites) recur_times <- c(0, recur_times)
  
  return(recur_times)
}

# simulate sequence of primary infection and bite times
simulate_patient <- function(LAMBDA, TIME_STEP, N_AGE, N_OBS, NU, R, ETA, PREL) {
  bite_times <- simulate_bite_times(rev(rev(LAMBDA)[1:(N_AGE+N_OBS)]), TIME_STEP, N_AGE+N_OBS)
  
  lapply(bite_times, function(x) x + simulate_recurrence_times(NU, R, ETA, PREL)) %>% 
    do.call(c, .) %>% sort %>% return
}

# N_AGES: vector of individual patient ages to simulate
simulate_cohort <- function(LAMBDA, TIME_STEP, N_AGES, N_OBS, NU, R, ETA, PREL, U, SHIFT_AGE) {
  t(sapply(N_AGES, function(N_AGE) {
    LAMBDA_ADJUSTED <- c(head(LAMBDA, -(N_OBS+N_AGE-SHIFT_AGE))*U, tail(LAMBDA, N_OBS+N_AGE-SHIFT_AGE))
    inf <- simulate_patient(LAMBDA_ADJUSTED, TIME_STEP, N_AGE, N_OBS, NU, R, ETA, PREL) - N_AGE*TIME_STEP
    inf <- (inf[inf>=0 & inf<=N_OBS*TIME_STEP])%/%TIME_STEP+1
    y <- rep(0, N_OBS); y[inf] <- 1; return(y)}))
}

# simulate masking --- either due to prophylaxis or immunity
# PCLIN should be a vector (masking probability by individual/age)
simulate_masking <- function(inf_matrix, PCLIN, PROPHYLAXIS_WINDOW, 
                             BUNCHING_WINDOW, MASKING_MATRIX) {
  inf_matrix <- inf_matrix * MASKING_MATRIX
  
  for (i in 1:nrow(inf_matrix)) {
    inf_times <- which(inf_matrix[i,]==1)
    masked_times <- which(is.na(inf_matrix[i,]))
    nonclinical <- inf_times[runif(length(inf_times))>PCLIN[i]]
    inf_matrix[i, nonclinical] <- 0
    
    if (PROPHYLAXIS_WINDOW>0) {
      for (j in 1:ncol(inf_matrix)) {
        if (!is.na(inf_matrix[i, j]) & inf_matrix[i, j]==1) {
          inf_matrix[i,  subset(j+(1:PROPHYLAXIS_WINDOW), j+(1:PROPHYLAXIS_WINDOW)<=ncol(inf_matrix))] <- NA
          
          check <- any(inf_matrix[i,  subset(j+PROPHYLAXIS_WINDOW+(1:BUNCHING_WINDOW), 
                                             j+PROPHYLAXIS_WINDOW+(1:BUNCHING_WINDOW)<=ncol(inf_matrix))]==1)
          check <- ifelse(is.na(check), FALSE, check)
          if (check) {
            inf_matrix[i,  subset(j+PROPHYLAXIS_WINDOW+BUNCHING_WINDOW, 
                                  j+PROPHYLAXIS_WINDOW+BUNCHING_WINDOW<=ncol(inf_matrix))] <- 1
            inf_matrix[i,  subset(j+PROPHYLAXIS_WINDOW+(1:(BUNCHING_WINDOW-1)), 
                                  j+PROPHYLAXIS_WINDOW+(1:(BUNCHING_WINDOW-1))<=ncol(inf_matrix))] <- 0
          }
        }
      }
    }
    inf_matrix[i, masked_times] <- NA 
  }
  return(inf_matrix)
}

# simulate a complete dataset with masking an an upper bound on the number of 
# infections experienced by each child
simulate_dataset <- function(LAMBDA, TIME_STEP, N_AGES, N_OBS, NU, R, ETA, PREL, 
                             PROPHYLAXIS_WINDOW, BUNCHING_WINDOW, PCLIN, MASKING_MATRIX, U, SHIFT_AGE) {
  
  all_inf <- simulate_cohort(LAMBDA, TIME_STEP, N_AGES, N_OBS, NU, R, ETA, PREL, U, SHIFT_AGE)
  inf_matrix <- simulate_masking(all_inf, PCLIN, PROPHYLAXIS_WINDOW, BUNCHING_WINDOW, MASKING_MATRIX)
  
  return(inf_matrix)
}