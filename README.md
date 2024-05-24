This repository contains code and a mininum analysis dataset to accompany the manuscript **Modelling the within-host dynamics of *Plasmodium vivax*
hypnozoite activation: an analysis of the SPf66 vaccine trial** by Mehra et al, 2024. This is available as a preprint on medRxiv: <https://www.medrxiv.org/content/10.1101/2024.01.03.24300707v1>.

The analysis can be reproduced by running the R Markdown files as follows:
* *Spf66_data_processing.Rmd*: this generates summary metrics for the Spf66 trial; performs adjustments for post-treatment prophylaxis and treatment failure; and pre-processes data for model calibration.
* *Spf66_model_calibration.Rmd*: this fits a mathematical model (predicated on the exponential clock model of hypnozoite activation) to the temporal sequence of symptomatic vivax episodes in the Spf66 vaccine trial; these are the results presented in the main text of the manuscript.
* *Spf66_supplementary_fits.Rmd*: this performs supplementary model fits to interrogate the consequences of certain forms of model misspecification, as addressed in the appendix.

A minimum analysis dataset, including patient metadata and a summary of malaria consultation records, can be found in the directory *Spf66_data_cleaned*.

The following packages are required to run the RMarkdown files:
* *Handling data*: dplyr, tidyr, readr, lubridate, IRanges
* *Model fitting*: Rcpp, RcppArmadillo, parallel
* *Statistical tests*: survival, MatchIt, PMCMRplus
* *Visualisation*: ggplot2, cowplot, ggfortify, GGally, knitr, kableExtra

Model fitting is relatively computationally expensive. After parallelisation, fitting the base model takes around 8 to 10 hours. Collectively, the simulation studies and supplementary fits can take around a week to run on a standard laptop. For completeness, all posteriors have been saved in .rds format as intermediate files. The analysis can be re-run from scratch by setting the variable RUN_ALL to TRUE in each RMarkdown file. If the variable RUN_ALL is set to FALSE, intermediate files are used to generate model summaries and perform posterior predictive checks. 