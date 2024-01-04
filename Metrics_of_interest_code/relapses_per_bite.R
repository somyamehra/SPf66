detectable_relapses_per_bite <- function(NU, PREL, ETA, T_MASK) {
  PROB_MASKED <- exp(-T_MASK*ETA)
  
  # PGFs for number of activatable hypnozoites remaining after ith relapse window
  G0 <- function(z) {1/(1+NU*PREL*(1-z))}
  G1 <- function(z) {G0(0) + (G0(1+(z-1)*PROB_MASKED) - G0(0))/(1+(z-1)*PROB_MASKED)}
  G2 <- function(z) {G1(0) + (G1(1+(z-1)*PROB_MASKED) - G1(0))/(1+(z-1)*PROB_MASKED)}
  G3 <- function(z) {G2(0) + (G2(1+(z-1)*PROB_MASKED) - G2(0))/(1+(z-1)*PROB_MASKED)}
  G4 <- function(z) {G3(0) + (G3(1+(z-1)*PROB_MASKED) - G3(0))/(1+(z-1)*PROB_MASKED)}
  G5 <- function(z) {G4(0) + (G4(1+(z-1)*PROB_MASKED) - G4(0))/(1+(z-1)*PROB_MASKED)}
  G6 <- function(z) {G5(0) + (G5(1+(z-1)*PROB_MASKED) - G5(0))/(1+(z-1)*PROB_MASKED)}
  G7 <- function(z) {G6(0) + (G6(1+(z-1)*PROB_MASKED) - G6(0))/(1+(z-1)*PROB_MASKED)}
  G8 <- function(z) {G7(0) + (G7(1+(z-1)*PROB_MASKED) - G7(0))/(1+(z-1)*PROB_MASKED)}
  G9 <- function(z) {G8(0) + (G8(1+(z-1)*PROB_MASKED) - G8(0))/(1+(z-1)*PROB_MASKED)}
  G10 <- function(z) {G9(0) + (G9(1+(z-1)*PROB_MASKED) - G9(0))/(1+(z-1)*PROB_MASKED)}
  G11 <- function(z) {G10(0) + (G10(1+(z-1)*PROB_MASKED) - G10(0))/(1+(z-1)*PROB_MASKED)}
  G12 <- function(z) {G11(0) + (G11(1+(z-1)*PROB_MASKED) - G11(0))/(1+(z-1)*PROB_MASKED)}
  G13 <- function(z) {G12(0) + (G12(1+(z-1)*PROB_MASKED) - G12(0))/(1+(z-1)*PROB_MASKED)}
  G14 <- function(z) {G13(0) + (G13(1+(z-1)*PROB_MASKED) - G13(0))/(1+(z-1)*PROB_MASKED)}
  G15 <- function(z) {G14(0) + (G14(1+(z-1)*PROB_MASKED) - G14(0))/(1+(z-1)*PROB_MASKED)}
  
  # number of relapses
  detected_relapses_pmf <- 
    1-c(G0(0), G1(0), G2(0), G3(0), G4(0), G5(0), G6(0), G7(0), G8(0), G9(0), G10(0), 
        G11(0), G12(0), G13(0), G14(0), G15(0))
  
  all_relapses_pmf <- 1-cumsum(sapply(0:15, function(x) {(NU*PREL)^x/(1+NU*PREL)^(1+x)}))
  
  return(data.frame(n_relapses=1:16, activated=all_relapses_pmf, 
                    detected=detected_relapses_pmf))
  
} 

inter_relapse_cdfs <- function(NU, PREL, ETA, T_MASK, T_VALS) {
  PROB_MASKED <- exp(-T_MASK*ETA)
  EXP_T_VALS <- exp(-ETA*T_VALS)
  
  F1 <- function(z1) {1/(1+NU*PREL*(1-z1))}
  
  F2 <- function(z1, z2) {
    F1(0) + (F1(z1*(1+(z2-1)*PROB_MASKED)) - F1(0))/(1+(z2-1)*PROB_MASKED)}
  
  F3 <- function(z1, z2, z3) {
    F2(z1, 0) + (F2(z1, z2*(1+(z3-1)*PROB_MASKED)) - F2(z1, 0))/(1+(z3-1)*PROB_MASKED)}
  
  F4 <- function(z1, z2, z3, z4) {
    F3(z1, z2, 0) + (F3(z1, z2, z3*(1+(z4-1)*PROB_MASKED)) - F3(z1, z2, 0))/(1+(z4-1)*PROB_MASKED)}
  
  F5 <- function(z1, z2, z3, z4, z5) {
    F4(z1, z2, z3, 0) + (F4(z1, z2, z3, z4*(1+(z5-1)*PROB_MASKED)) - F4(z1, z2, z3, 0))/(1+(z5-1)*PROB_MASKED)}
  
  F6 <- function(z1, z2, z3, z4, z5, z6) {
    F5(z1, z2, z3, z4, 0) + (F5(z1, z2, z3, z4, z5*(1+(z6-1)*PROB_MASKED)) - F5(z1, z2, z3, z4, 0))/(1+(z6-1)*PROB_MASKED)}
  
  F7 <- function(z1, z2, z3, z4, z5, z6, z7) {
    F6(z1, z2, z3, z4, z5, 0) + (F6(z1, z2, z3, z4, z5, z6*(1+(z7-1)*PROB_MASKED)) - F6(z1, z2, z3, z4, z5, 0))/(1+(z7-1)*PROB_MASKED)}
  
  F8 <- function(z1, z2, z3, z4, z5, z6, z7, z8) {
    F7(z1, z2, z3, z4, z5, z6, 0) + (F7(z1, z2, z3, z4, z5, z6, z7*(1+(z8-1)*PROB_MASKED)) - F7(z1, z2, z3, z4, z5, z6, 0))/(1+(z8-1)*PROB_MASKED)}
  
  F9 <- function(z1, z2, z3, z4, z5, z6, z7, z8, z9) {
    F8(z1, z2, z3, z4, z5, z6, z7, 0) + (F8(z1, z2, z3, z4, z5, z6, z7, z8*(1+(z9-1)*PROB_MASKED)) - F8(z1, z2, z3, z4, z5, z6, z7, 0))/(1+(z9-1)*PROB_MASKED)}
  
  
  inter_relapse_cdf <- list()
  
  inter_relapse_cdf[["2 detectable relapses"]] <- 
    data.frame(time=T_VALS,
               relapse_1=(F3(EXP_T_VALS, 1, 0)-F3(EXP_T_VALS, 0, 0))/(F3(1, 1, 0)-F3(1, 0, 0)),
               relapse_2=(F3(1, EXP_T_VALS, 0)-F3(1, 0, 0))/(F3(1, 1, 0)-F3(1, 0, 0)))
  
  
  inter_relapse_cdf[["4 detectable relapses"]] <- 
    data.frame(time=T_VALS,
               relapse_1=(F5(EXP_T_VALS, 1, 1, 1, 0)-F5(EXP_T_VALS, 1, 1, 0, 0))/(F5(1, 1, 1, 1, 0)-F5(1, 1, 1, 0, 0)),
               relapse_2=(F5(1, EXP_T_VALS, 1, 1, 0)-F5(1, EXP_T_VALS, 1, 0, 0))/(F5(1, 1, 1, 1, 0)-F5(1, 1, 1, 0, 0)),
               relapse_3=(F5(1, 1, EXP_T_VALS, 1, 0)-F5(1, 1, EXP_T_VALS, 0, 0))/(F5(1, 1, 1, 1, 0)-F5(1, 1, 1, 0, 0)),
               relapse_4=(F5(1, 1, 1, EXP_T_VALS, 0)-F5(1, 1, 1, 0, 0))/(F5(1, 1, 1, 1, 0)-F5(1, 1, 1, 0, 0)))
  
  inter_relapse_cdf[["6 detectable relapses"]] <- 
    x <- data.frame(time=T_VALS,
                    relapse_1=(F7(EXP_T_VALS, 1, 1, 1, 1, 1, 0)-F7(EXP_T_VALS, 1, 1, 1, 1, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)),
                    relapse_2=(F7(1, EXP_T_VALS, 1, 1, 1, 1, 0)-F7(1, EXP_T_VALS, 1, 1, 1, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)),
                    relapse_3=(F7(1, 1, EXP_T_VALS, 1, 1, 1, 0)-F7(1, 1, EXP_T_VALS, 1, 1, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)),
                    relapse_4=(F7(1, 1, 1, EXP_T_VALS, 1, 1, 0)-F7(1, 1, 1, EXP_T_VALS, 1, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)),
                    relapse_5=(F7(1, 1, 1, 1, EXP_T_VALS, 1, 0)-F7(1, 1, 1, 1, EXP_T_VALS, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)),
                    relapse_6=(F7(1, 1, 1, 1, 1, EXP_T_VALS, 0)-F7(1, 1, 1, 1, 1, 0, 0))/
                      (F7(1, 1, 1, 1, 1, 1, 0)-F7(1, 1, 1, 1, 1, 0, 0)))
  
  inter_relapse_cdf[["8 detectable relapses"]] <- 
    x <- data.frame(time=T_VALS,
                    relapse_1=(F9(EXP_T_VALS, 1, 1, 1, 1, 1, 1, 1, 0)-F9(EXP_T_VALS, 1, 1, 1, 1, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_2=(F9(1, EXP_T_VALS, 1, 1, 1, 1, 1, 1, 0)-F9(1, EXP_T_VALS, 1, 1, 1, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_3=(F9(1, 1, EXP_T_VALS, 1, 1, 1, 1, 1, 0)-F9(1, 1, EXP_T_VALS, 1, 1, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_4=(F9(1, 1, 1, EXP_T_VALS, 1, 1, 1, 1, 0)-F9(1, 1, 1, EXP_T_VALS, 1, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_5=(F9(1, 1, 1, 1, EXP_T_VALS, 1, 1, 1, 0)-F9(1, 1, 1, 1, EXP_T_VALS, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_6=(F9(1, 1, 1, 1, 1, EXP_T_VALS, 1, 1, 0)-F9(1, 1, 1, 1, 1, EXP_T_VALS, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_7=(F9(1, 1, 1, 1, 1, 1, EXP_T_VALS, 1, 0)-F9(1, 1, 1, 1, 1, 1, EXP_T_VALS, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)),
                    relapse_8=(F9(1, 1, 1, 1, 1, 1, 1, EXP_T_VALS, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0))/
                      (F9(1, 1, 1, 1, 1, 1, 1, 1, 0)-F9(1, 1, 1, 1, 1, 1, 1, 0, 0)))
  
  
  inter_relapse_cdf <- lapply(inter_relapse_cdf, function(x) {
    x %>% reshape2::melt(id="time") %>% 
      mutate(time=time+T_MASK*(variable!="relapse_1")) %>%
      bind_rows(data.frame(time=0, variable=unique(.$variable), value=1))}) %>%
    bind_rows(.id="n_detected")
}