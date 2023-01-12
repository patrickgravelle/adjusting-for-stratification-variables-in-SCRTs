# Adjusting For Stratification Variables in Cluster Randomized Trials

## About

Research article by Patrick Gravelle and Dr. Roee Gutman completed as part of Patrick Gravelle's doctoral dissertation for Brown's Biostatistics Department. Article submitted to Statistics in Medicine journal.

## Abstract

Randomized controlled trials (RCTs) are the gold-standard for estimating the effects of interventions. While adjustment for covariates in RCTs leads to more powerful statistical procedures, multiple systematic reviews have shown that many RCTs do not adjust for covariates beyond the intervention. In cluster randomized trials (CRTs) groups of subjects, rather than individuals, are randomly assigned to interventions. Regression adjustments in CRTs rely on more complex modelling because they commonly comprise of individual-level and cluster-level covariates as well as adjustments for the correlations of individuals within groups. Multilevel generalized linear models (MGLM) and generalized estimating equations (GEE) have been proposed as possible analysis methods for individual-level data in CRTs. There is limited literature on the exact specification of these regression models with CRTs and only limited comparisons between them. Moreover, there is limited evidence on the operating characteristics of GEE and MGLM in CRTs with stratified randomization. Using extensive simulations with binary and continuous outcomes, we compare the performance of different specifications of MGLM and GEE models. Among the models considered, the model proposed by Neuhaus and Kalbfleisch (1998) is a valid statistical procedure that is most precise and accurate for continuous and binary outcomes. GEE procedures do not offer improvements compared to MGLMs with worse performance for binary outcomes.

## Code

See the accompanying R code to replicate the simulations. Additional simulations can be constructed using the configurations outlined in the paper along with the code provided.
