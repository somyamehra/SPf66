req_packages <- c("dplyr", "tidyr", "readr", "lubridate", 
                  "BiocManager", "Rcpp", "RcppArmadillo", "parallel", 
                  "survival", "MatchIt", "PMCMRplus",
                  "ggplot2", "cowplot", "ggfortify", "GGally", "knitr", "kableExtra")

sapply(req_packages, function(x) {
  if (!(x %in% installed.packages()[, "Package"])) install.packages(x, dependencies = TRUE)
  require(x, character.only = TRUE)
})

if (!("IRanges" %in% installed.packages()[, "Package"])) BiocManager::install("IRanges")