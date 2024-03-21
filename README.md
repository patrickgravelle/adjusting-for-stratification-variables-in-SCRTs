# Adjusting For Stratification Variables in Cluster Randomized Trials

## About

Research article by Patrick Gravelle and Dr. Roee Gutman completed as part of Patrick Gravelle's doctoral dissertation for Brown's Biostatistics Department. Article submitted for publication. This project was supported by the IMbedded Pragmatic Alzheimer’s Disease (AD) and AD-Related Dementias (AD/ADRD) Clinical Trials (IMPACT) Collaboratory Design and Statistics Core which is NIA grant U54AG063546.

## Abstract

Randomized controlled trials (RCTs) are the gold-standard for estimating the effects of interventions. While adjustment for covariates in RCTs leads to more powerful statistical procedures, multiple systematic reviews have shown that many RCTs do not adjust for covariates beyond the intervention. In cluster randomized trials (CRTs) groups of subjects, rather than individuals, are randomly assigned to interventions. When there exists a variable believed to be prognostic of the outcome of interest, stratifying the randomization by this variable may lead to improved balance between intervention arms. In the analysis phase of a stratified CRT (SCRT), generalized linear multilevel models (GLMM) and generalized estimating equations (GEE) have been proposed to adjust for individual-level and cluster-level covariates as well as for the correlations between individuals within the same group. When a SCRT is implemented, the exact specification of GLMM and GEE models that account for stratification variables have not been extensively examined. We describe possible specifications of GLMMs and GEE models for SCRTs. Using extensive simulations with continuous and binary outcomes, we compare the performance of these models. Based on the simulations, we identify four main conclusions. First, models that adjust for individual and cluster levels covariates have better precision than models that do not. Second, among the models considered, a model that adjusts for the within and between cluster association of baseline covariates is a valid statistical procedure that is most precise and accurate for continuous and binary outcomes. Third, GEE models should use bias-adjusted sampling variance to obtain valid statistical inference. Fourth, GEE models generally result in wider interval estimates compared to corresponding GLMMs with similar mean structure. 

## Code

See the accompanying R code to replicate the simulations. Additional simulations can be constructed using the configurations outlined in the paper along with the code provided.
