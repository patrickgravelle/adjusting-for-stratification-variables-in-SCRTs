
library(dplyr)
library(geepack)
library(rstan)
library(brms)
library(brmsmargins)
library(marginaleffects)

## Coverage, MSE, and Bias matrix
crt_design4_missp_RE_Gamma <- matrix(NA, nrow = 243, ncol = 28)

# counter4later_matrix

## Indices for placing values properly
index <- 0

n_units <- 20

n_sims <- 100




# 2 strata 

## Define the number of clusters, units, and simulations
n_clusters <- rep(c(48,72,96), 27)

## Define the effect of the treatment
treat_effect <- rep(rep(c(0.25,1,0.5), each = 9), 3)

## Define the effect of the proportion of race
prop_effect <- rep(c(1,1.5,2), each = 27)

cluster_var <- rep(c(rep(0.5, 3), rep(2, 3), rep(4, 3)), 9)

cont_effect <- 10


seeds <- (1:(length(cluster_var)*n_sims*3) + (length(cluster_var)*n_sims*3)*20)



priorset <- c(set_prior("normal(0,5)", class = "b"),
              set_prior("student_t(3,0,5)", class = "sd"))


for (j in 1:length(cluster_var)){
  
  ## Matrix to hold the treatment effects and standard errors
  treatment_diff <- matrix(NA, nrow = n_sims, ncol = 15)
  
  ## indexing through the coverage matrix
  index <- index + 1
  
  singular_check <- matrix(NA, nrow = n_sims, ncol = 5)
  
  
  ################################
  
  # Iterate through the simulations
  
  for (m in 1:n_sims){
    
    set.seed(seeds[((j-1)*3*n_sims + m + 0*n_sims)])
    
    # random mean for each cluster
    cluster_mean <- rep(rgamma(n_clusters[j], cluster_var[j], 1), each = n_units)
    
    # creating the cluster vector
    cluster <- rep(1:n_clusters[j], each = n_units)
    
    # extra_covariate <- rep(-1:(n_clusters[j]-2), each = n_units)
    
    # random error for each observation
    noise <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    # dataframe for the factors
    full_data <- data.frame(cluster_mean, cluster,   noise)
    
    names(full_data) <- c("cluster_mean", "cluster",   "noise")
    
    # proportion race variable for each cluster
    prop <- runif(n_clusters[j], 0, 1)
    
    # full length proportion vector (same value for people in the same cluster)
    full_data$proportion <- rep(prop, each = n_units)
    
    # order the dataset by proportion variable
    full_data <- full_data[order(full_data$proportion), ]
    
    # individual race (0 or 1) using the cluster proportion as the probability of being race
    full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
    
    
    # add strata identifier vector
    prop_strat <- c(1,2)
    full_data$prop_strata <- rep(prop_strat, each = ((n_clusters[j]/2)*n_units))
    
    ## Continuous predictor
    full_data$continuous_pred <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
    
    # order the dataset by cluster mean variable
    full_data <- full_data[order(full_data$cont_ave), ]
    
    
    # # add continuous strata identifier vector
    cont_strat <- c(1,2)
    full_data$cont_strata <- rep(cont_strat, each = ((n_clusters[j]/2)*n_units))
    
    
    ## Add strata vectors into the dataset
    # full_data <- data.frame(full_data, cont_strata)
    
    ## New Strata based on both cont and prop strata
    full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                   ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                          ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 4) ))
    
    
    clust1 <- nrow(full_data[full_data$all_strata==1, ]) / n_units
    clust2 <- nrow(full_data[full_data$all_strata==2, ]) / n_units
    clust3 <- nrow(full_data[full_data$all_strata==3, ]) / n_units
    clust4 <- nrow(full_data[full_data$all_strata==4, ]) / n_units
    
    
    ## Stratify data by proportion race (here we use two strata)
    ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
    ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
    treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
    treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
    treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
    treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
    
    ## Combine the treatment vectors
    # full treatment vector
    treated <- c(treated1, treated2, treated3, treated4)
    
    # order data by all strata
    full_data <- full_data[order(full_data$all_strata), ]
    
    # Replicate treatment assignment for the number of people per cluster
    full_data$treatment <- rep(treated, each = n_units)
    
    
    
    ## Generate our outcome by taking the sum of the following:
    # Cluster mean
    # treatment assignment (0 or 1) times the defined treatment effect
    # proportion race (cluster level) times the defined proportion effect
    # race (0 or 1)
    # noise
    
    linear_predictor <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment -
      (prop_effect[j])*full_data$proportion - 
      full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
    
    probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
    
    
    full_data$outcomes <- rbinom(n_clusters[j]*n_units, 1 , probability)
    
    ## For one of the literature models we create the "within cluster" correlation variable defined as follows
    ## individual value of race divided by the mean (here I used the proportion rather than the specific cluster mean)
    ## This would only effect this model's specific results
    full_data$within <- full_data$race - full_data$proportion
    
    full_data$within_cont <- full_data$continuous_pred - full_data$cont_ave
    
    # Now we define our models
    
    # Now we define our models
    if (j != 1 | m != 1){
      
      # Now we define our models
      small_model <- update(small_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      strat_model <- update(strat_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      large_model <- update(large_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      kalbfleisch <- update(kalbfleisch_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      teerenstra <- update(teerenstra_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
    } else {
      
      small_model_1 <- brm(outcomes ~ treatment + race + (1|cluster), 
                         data = full_data, family = "bernoulli", prior = priorset,
                         warmup = 500, iter = 3000, chains = 1,
                         control = list(adapt_delta = 0.9), seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      strat_model_1 <- brm(outcomes ~ treatment + race + all_strata + (1|cluster), 
                         data = full_data, family = "bernoulli", prior = priorset,
                         warmup = 500, iter = 3000, chains = 1,
                         control = list(adapt_delta = 0.9), seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      large_model_1 <- brm(outcomes ~ treatment + race + proportion + exp(continuous_pred) + 
                           cont_ave + (1|cluster), 
                         data = full_data, family = "bernoulli", prior = priorset,
                         warmup = 500, iter = 3000, chains = 1,
                         control = list(adapt_delta = 0.9), seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      kalbfleisch_1 <- brm(outcomes ~ treatment + proportion + within + cont_ave + within_cont + (1|cluster), 
                         data = full_data, family = "bernoulli", prior = priorset,
                         warmup = 500, iter = 3000, chains = 1,
                         control = list(adapt_delta = 0.9), seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      teerenstra_1 <- brm(outcomes ~ treatment + (1|cluster), 
                        data = full_data, family = "bernoulli", prior = priorset,
                        warmup = 500, iter = 3000, chains = 1,
                        control = list(adapt_delta = 0.9), seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)])
      
      
      small_model <- small_model_1
      strat_model <- strat_model_1
      large_model <- large_model_1
      kalbfleisch <- kalbfleisch_1
      teerenstra <- teerenstra_1
    }
    
    
    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion + exp(continuous_pred) + cont_ave, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    singular_check[m,1] <- is.null(small_model)
    singular_check[m,2] <- is.null(strat_model)
    singular_check[m,3] <- is.null(large_model)
    singular_check[m,4] <- is.null(kalbfleisch)
    singular_check[m,5] <- is.null(teerenstra)
    
    if ((sum(singular_check[m,]) == 0) & (sum(is.null(gee_true) + is.null(gee_under)) == 0)){
      
      # small model
      marginal_model2 <- brmsmargins(
        small_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)], CI = 0.95)
      
      treatment_diff[m,2] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,3] <- sd(marginal_model2$Contrasts)
      
      # Repeat the same steps above for every model
      
      # true model
      marginal_model2 <- brmsmargins(
        large_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)], CI = 0.95)
      
      treatment_diff[m,4] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,5] <- sd(marginal_model2$Contrasts)
      
      # strat model
      marginal_model2 <- brmsmargins(
        strat_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)], CI = 0.95)
      
      treatment_diff[m,6] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,7] <- sd(marginal_model2$Contrasts)
      
      
      # kalbfleisch model
      marginal_model2 <- brmsmargins(
        kalbfleisch,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)], CI = 0.95)
      
      treatment_diff[m,8] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,9] <- sd(marginal_model2$Contrasts)
      
      
      # teerenstra model 
      marginal_model2 <- brmsmargins(
        teerenstra,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 0*n_sims)], CI = 0.95)
      
      treatment_diff[m,10] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,11] <- sd(marginal_model2$Contrasts)
      
      
      # gee true model
      marginal_model2 <- marginaleffects(gee_true)
      
      treatment_diff[m,12] <- summary(marginal_model2)[1,4]
      treatment_diff[m,13] <- summary(marginal_model2)[1,5]
      
      # gee under model
      marginal_model2 <- marginaleffects(gee_under)
      
      treatment_diff[m,14] <- summary(marginal_model2)[1,4]
      treatment_diff[m,15] <- summary(marginal_model2)[1,5]
      
    }
  }
  
  # random mean for each cluster
  cluster_mean <- rep(rgamma(18000, cluster_var[j], 1), each = 100)
  
  # creating the cluster vector
  cluster <- rep(1:18000, each = 100)
  
  # extra_covariate <- rep(-1:(18000-2), each = 100)
  
  # random error for each observation
  noise <- rnorm(18000*100, 0, 1)
  
  # dataframe for the factors
  full_data <- data.frame(cluster_mean, cluster,   noise)
  
  names(full_data) <- c("cluster_mean", "cluster",   "noise")
  
  # proportion race variable for each cluster
  prop <- runif(18000, 0, 1)
  
  # full length proportion vector (same value for people in the same cluster)
  full_data$proportion <- rep(prop, each = 100)
  
  # order the dataset by proportion variable
  full_data <- full_data[order(full_data$proportion), ]
  
  # individual race (0 or 1) using the cluster proportion as the probability of being race
  full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
  
  
  # add strata identifier vector
  prop_strat <- c(1,2)
  full_data$prop_strata <- rep(prop_strat, each = ((18000/2)*100))
  
  ## Continuous predictor
  full_data$continuous_pred <- rnorm(18000*100, 0, 1)
  
  full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
  
  # order the dataset by cluster mean variable
  full_data <- full_data[order(full_data$cont_ave), ]
  
  
  # # add continuous strata identifier vector
  cont_strat <- c(1,2)
  full_data$cont_strata <- rep(cont_strat, each = ((18000/2)*100))
  
  
  ## Add strata vectors into the dataset
  # full_data <- data.frame(full_data, cont_strata)
  
  ## New Strata based on both cont and prop strata
  full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                 ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                        ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 4) ))
  
  
  clust1 <- nrow(full_data[full_data$all_strata==1, ]) / 100
  clust2 <- nrow(full_data[full_data$all_strata==2, ]) / 100
  clust3 <- nrow(full_data[full_data$all_strata==3, ]) / 100
  clust4 <- nrow(full_data[full_data$all_strata==4, ]) / 100
  
  
  ## Stratify data by proportion race (here we use two strata)
  ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
  ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
  treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
  treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
  treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
  treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
  
  ## Combine the treatment vectors
  # full treatment vector
  treated <- c(treated1, treated2, treated3, treated4)
  
  # order data by all strata
  full_data <- full_data[order(full_data$all_strata), ]
  
  # Replicate treatment assignment for the number of people per cluster
  full_data$treatment <- rep(treated, each = 100)
  
  
  
  ## Generate our outcome by taking the sum of the following:
  # Cluster mean
  # treatment assignment (0 or 1) times the defined treatment effect
  # proportion race (cluster level) times the defined proportion effect
  # race (0 or 1)
  # noise
  
  linear_predictor <- full_data$cluster_mean +
    (treat_effect[j])*full_data$treatment -
    (prop_effect[j])*full_data$proportion - 
    full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
  
  probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
  
  
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  treatment_diff[,1] <- rep((mean(probability[which(full_data$treatment == 1)]) - mean(probability[which(full_data$treatment == 0)])), n_sims)
  
  # counter4later_total
  
  ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  treatment_diff <- na.omit(treatment_diff)
  
  ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  crt_design4_missp_RE_Gamma[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
                                              treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
                                              treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
                                              treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
                                              treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
                                              treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                              treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                              treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  crt_design4_missp_RE_Gamma[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  crt_design4_missp_RE_Gamma[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  
  crt_design4_missp_RE_Gamma[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  
  crt_design4_missp_RE_Gamma[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  
  crt_design4_missp_RE_Gamma[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  
  crt_design4_missp_RE_Gamma[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_design4_missp_RE_Gamma[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_design4_missp_RE_Gamma[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  crt_design4_missp_RE_Gamma[index, 22] <- mean(2*1.96*treatment_diff[,3])
  
  crt_design4_missp_RE_Gamma[index, 23] <- mean(2*1.96*treatment_diff[,7])
  
  crt_design4_missp_RE_Gamma[index, 24] <- mean(2*1.96*treatment_diff[,5])
  
  crt_design4_missp_RE_Gamma[index, 25] <- mean(2*1.96*treatment_diff[,9])
  
  crt_design4_missp_RE_Gamma[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_design4_missp_RE_Gamma[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_design4_missp_RE_Gamma[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  
  
  
  
  
  
  # 3 Strata
  
  
  treatment_diff <- matrix(NA, nrow = n_sims, ncol = 15)
  
  index <- index + 1
  
  singular_check <- matrix(NA, nrow = n_sims, ncol = 5)
  
  
  ################################
  
  for (m in 1:n_sims){
    
    set.seed(seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    # predictors for the model
    cluster_mean <- rep(rgamma(n_clusters[j], cluster_var[j], 1), each = n_units)
    
    cluster <- rep(1:n_clusters[j], each = n_units)
    
    # extra_covariate <- rep(-1:(n_clusters[j]-2), each = n_units)
    
    noise <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    full_data <- data.frame(cluster_mean, cluster,   noise)
    
    names(full_data) <- c("cluster_mean", "cluster",   "noise")
    
    # proportion variable on which to stratify data
    prop <- runif(n_clusters[j], 0, 1)
    
    # order the dataset by proportion variable
    full_data$proportion <- rep(prop, each = n_units)
    full_data <- full_data[order(full_data$proportion), ]
    
    full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
    
    # add strata identifier vector
    prop_strat <- c(1,2,3)
    prop_strata <- rep(prop_strat, each = ((n_clusters[j]/3)*n_units))
    
    ## Add strata vectors into the dataset
    full_data <- data.frame(full_data, prop_strata)
    
    ## Continuous predictor
    full_data$continuous_pred <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
    
    # order the dataset by cluster mean variable
    full_data <- full_data[order(full_data$cont_ave), ]
    
    
    # # add continuous strata identifier vector
    cont_strat <- c(1,2)
    full_data$cont_strata <- rep(cont_strat, each = ((n_clusters[j]/2)*n_units))
    
    
    ## New Strata based on both cont and prop strata
    full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                   ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                          ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 
                                                 ifelse((full_data$prop_strata==2 & full_data$cont_strata==2), 4, 
                                                        ifelse((full_data$prop_strata==3 & full_data$cont_strata==1), 5, 6))) ))
    
    
    clust1 <- nrow(full_data[full_data$all_strata==1, ]) / n_units
    clust2 <- nrow(full_data[full_data$all_strata==2, ]) / n_units
    clust3 <- nrow(full_data[full_data$all_strata==3, ]) / n_units
    clust4 <- nrow(full_data[full_data$all_strata==4, ]) / n_units
    clust5 <- nrow(full_data[full_data$all_strata==5, ]) / n_units
    clust6 <- nrow(full_data[full_data$all_strata==6, ]) / n_units
    
    
    
    ## Stratify data by proportion race (here we use two strata)
    ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
    ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
    treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
    treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
    treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
    treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
    treated5 <- sample(c(rep(0, round(clust5/2)), rep(1, clust5 - round(clust5/2))), replace = F)
    treated6 <- sample(c(rep(0, round(clust6/2)), rep(1, clust6 - round(clust6/2))), replace = F)
    
    ## Combine the treatment vectors
    # full treatment vector
    treated <- c(treated1, treated2, treated3, treated4, treated5, treated6)
    
    # order data by all strata
    full_data <- full_data[order(full_data$all_strata), ]
    
    # Replicate treatment assignment for the number of people per cluster
    full_data$treatment <- rep(treated, each = n_units)
    
    linear_predictor <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment -
      (prop_effect[j])*full_data$proportion - 
      full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
    
    probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
    
    full_data$outcomes <- rbinom(n_clusters[j]*n_units, 1 , probability)
    
    full_data$within <- full_data$race - full_data$proportion
    
    full_data$within_cont <- full_data$continuous_pred - full_data$cont_ave
    
    
    # Now we define our models
    small_model <- update(small_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    strat_model <- update(strat_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    large_model <- update(large_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    kalbfleisch <- update(kalbfleisch_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    teerenstra <- update(teerenstra_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    
    
    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion + exp(continuous_pred) + cont_ave, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    singular_check[m,1] <- is.null(small_model)
    singular_check[m,2] <- is.null(strat_model)
    singular_check[m,3] <- is.null(large_model)
    singular_check[m,4] <- is.null(kalbfleisch)
    singular_check[m,5] <- is.null(teerenstra)
    
    if ((sum(singular_check[m,]) == 0) & (sum(is.null(gee_true) + is.null(gee_under)) == 0)){
      
      # small model
      marginal_model2 <- brmsmargins(
        small_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)], CI = 0.95)
      
      treatment_diff[m,2] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,3] <- sd(marginal_model2$Contrasts)
      
      # Repeat the same steps above for every model
      
      # true model
      marginal_model2 <- brmsmargins(
        large_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)], CI = 0.95)
      
      treatment_diff[m,4] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,5] <- sd(marginal_model2$Contrasts)
      
      # strat model
      marginal_model2 <- brmsmargins(
        strat_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)], CI = 0.95)
      
      treatment_diff[m,6] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,7] <- sd(marginal_model2$Contrasts)
      
      
      # kalbfleisch model
      marginal_model2 <- brmsmargins(
        kalbfleisch,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)], CI = 0.95)
      
      treatment_diff[m,8] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,9] <- sd(marginal_model2$Contrasts)
      
      
      # teerenstra model 
      marginal_model2 <- brmsmargins(
        teerenstra,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 1*n_sims)], CI = 0.95)
      
      treatment_diff[m,10] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,11] <- sd(marginal_model2$Contrasts)
      
      
      # gee true model
      marginal_model2 <- marginaleffects(gee_true)
      
      treatment_diff[m,12] <- summary(marginal_model2)[1,4]
      treatment_diff[m,13] <- summary(marginal_model2)[1,5]
      
      # gee under model
      marginal_model2 <- marginaleffects(gee_under)
      
      treatment_diff[m,14] <- summary(marginal_model2)[1,4]
      treatment_diff[m,15] <- summary(marginal_model2)[1,5]
      
    }
  }
  
  
  # predictors for the model
  cluster_mean <- rep(rgamma(18000, cluster_var[j], 1), each = 100)
  
  cluster <- rep(1:18000, each = 100)
  
  # extra_covariate <- rep(-1:(18000-2), each = 100)
  
  noise <- rnorm(18000*100, 0, 1)
  
  full_data <- data.frame(cluster_mean, cluster,   noise)
  
  names(full_data) <- c("cluster_mean", "cluster",   "noise")
  
  # proportion variable on which to stratify data
  prop <- runif(18000, 0, 1)
  
  # order the dataset by proportion variable
  full_data$proportion <- rep(prop, each = 100)
  full_data <- full_data[order(full_data$proportion), ]
  
  full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
  
  # add strata identifier vector
  prop_strat <- c(1,2,3)
  prop_strata <- rep(prop_strat, each = ((18000/3)*100))
  
  ## Add strata vectors into the dataset
  full_data <- data.frame(full_data, prop_strata)
  
  ## Continuous predictor
  full_data$continuous_pred <- rnorm(18000*100, 0, 1)
  
  full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
  
  # order the dataset by cluster mean variable
  full_data <- full_data[order(full_data$cont_ave), ]
  
  
  # # add continuous strata identifier vector
  cont_strat <- c(1,2)
  full_data$cont_strata <- rep(cont_strat, each = ((18000/2)*100))
  
  
  ## New Strata based on both cont and prop strata
  full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                 ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                        ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 
                                               ifelse((full_data$prop_strata==2 & full_data$cont_strata==2), 4, 
                                                      ifelse((full_data$prop_strata==3 & full_data$cont_strata==1), 5, 6))) ))
  
  
  clust1 <- nrow(full_data[full_data$all_strata==1, ]) / 100
  clust2 <- nrow(full_data[full_data$all_strata==2, ]) / 100
  clust3 <- nrow(full_data[full_data$all_strata==3, ]) / 100
  clust4 <- nrow(full_data[full_data$all_strata==4, ]) / 100
  clust5 <- nrow(full_data[full_data$all_strata==5, ]) / 100
  clust6 <- nrow(full_data[full_data$all_strata==6, ]) / 100
  
  
  
  ## Stratify data by proportion race (here we use two strata)
  ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
  ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
  treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
  treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
  treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
  treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
  treated5 <- sample(c(rep(0, round(clust5/2)), rep(1, clust5 - round(clust5/2))), replace = F)
  treated6 <- sample(c(rep(0, round(clust6/2)), rep(1, clust6 - round(clust6/2))), replace = F)
  
  ## Combine the treatment vectors
  # full treatment vector
  treated <- c(treated1, treated2, treated3, treated4, treated5, treated6)
  
  # order data by all strata
  full_data <- full_data[order(full_data$all_strata), ]
  
  # Replicate treatment assignment for the number of people per cluster
  full_data$treatment <- rep(treated, each = 100)
  
  linear_predictor <- full_data$cluster_mean +
    (treat_effect[j])*full_data$treatment -
    (prop_effect[j])*full_data$proportion - 
    full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
  
  probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
  
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  treatment_diff[,1] <- rep((mean(probability[which(full_data$treatment == 1)]) - mean(probability[which(full_data$treatment == 0)])), n_sims)
  
  # counter4later_total
  
  ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  treatment_diff <- na.omit(treatment_diff)
  
  ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  crt_design4_missp_RE_Gamma[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
                                              treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
                                              treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
                                              treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
                                              treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
                                              treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                              treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                              treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  crt_design4_missp_RE_Gamma[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  crt_design4_missp_RE_Gamma[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  
  crt_design4_missp_RE_Gamma[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  
  crt_design4_missp_RE_Gamma[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  
  crt_design4_missp_RE_Gamma[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  
  crt_design4_missp_RE_Gamma[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_design4_missp_RE_Gamma[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_design4_missp_RE_Gamma[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  crt_design4_missp_RE_Gamma[index, 22] <- mean(2*1.96*treatment_diff[,3])
  
  crt_design4_missp_RE_Gamma[index, 23] <- mean(2*1.96*treatment_diff[,7])
  
  crt_design4_missp_RE_Gamma[index, 24] <- mean(2*1.96*treatment_diff[,5])
  
  crt_design4_missp_RE_Gamma[index, 25] <- mean(2*1.96*treatment_diff[,9])
  
  crt_design4_missp_RE_Gamma[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_design4_missp_RE_Gamma[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_design4_missp_RE_Gamma[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  
  
  # 4 Strata
  
  
  treatment_diff <- matrix(NA, nrow = n_sims, ncol = 15)
  
  index <- index + 1
  
  singular_check <- matrix(NA, nrow = n_sims, ncol = 5)
  
  
  ################################
  
  for (m in 1:n_sims){
    
    set.seed(seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    # predictors for the model
    cluster_mean <- rep(rgamma(n_clusters[j], cluster_var[j], 1), each = n_units)
    
    cluster <- rep(1:n_clusters[j], each = n_units)
    
    # extra_covariate <- rep(-1:(n_clusters[j]-2), each = n_units)
    
    noise <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    full_data <- data.frame(cluster_mean, cluster,   noise)
    
    names(full_data) <- c("cluster_mean", "cluster",   "noise")
    
    # proportion variable on which to stratify data
    prop <- runif(n_clusters[j], 0, 1)
    
    # order the dataset by proportion variable
    full_data$proportion <- rep(prop, each = n_units)
    full_data <- full_data[order(full_data$proportion), ]
    full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
    
    # add strata identifier vector
    prop_strat <- c(1,2,3,4)
    prop_strata <- rep(prop_strat, each = ((n_clusters[j]/4)*n_units))
    
    ## Add strata vectors into the dataset
    full_data <- data.frame(full_data, prop_strata)
    
    ## Continuous predictor
    full_data$continuous_pred <- rnorm(n_clusters[j]*n_units, 0, 1)
    
    full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
    
    # order the dataset by cluster mean variable
    full_data <- full_data[order(full_data$cont_ave), ]
    
    
    # # add continuous strata identifier vector
    cont_strat <- c(1,2)
    full_data$cont_strata <- rep(cont_strat, each = ((n_clusters[j]/2)*n_units))
    
    ## New Strata based on both cont and prop strata
    full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                   ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                          ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 
                                                 ifelse((full_data$prop_strata==2 & full_data$cont_strata==2), 4, 
                                                        ifelse((full_data$prop_strata==3 & full_data$cont_strata==1), 5, 
                                                               ifelse((full_data$prop_strata==3 & full_data$cont_strata==2), 6, 
                                                                      ifelse((full_data$prop_strata==4 & full_data$cont_strata==1), 7, 8))))) ))
    
    
    clust1 <- nrow(full_data[full_data$all_strata==1, ]) / n_units
    clust2 <- nrow(full_data[full_data$all_strata==2, ]) / n_units
    clust3 <- nrow(full_data[full_data$all_strata==3, ]) / n_units
    clust4 <- nrow(full_data[full_data$all_strata==4, ]) / n_units
    clust5 <- nrow(full_data[full_data$all_strata==5, ]) / n_units
    clust6 <- nrow(full_data[full_data$all_strata==6, ]) / n_units
    clust7 <- nrow(full_data[full_data$all_strata==7, ]) / n_units
    clust8 <- nrow(full_data[full_data$all_strata==8, ]) / n_units
    
    
    
    ## Stratify data by proportion race (here we use two strata)
    ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
    ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
    treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
    treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
    treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
    treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
    treated5 <- sample(c(rep(0, round(clust5/2)), rep(1, clust5 - round(clust5/2))), replace = F)
    treated6 <- sample(c(rep(0, round(clust6/2)), rep(1, clust6 - round(clust6/2))), replace = F)
    treated7 <- sample(c(rep(0, round(clust7/2)), rep(1, clust7 - round(clust7/2))), replace = F)
    treated8 <- sample(c(rep(0, round(clust8/2)), rep(1, clust8 - round(clust8/2))), replace = F)
    
    ## Combine the treatment vectors
    # full treatment vector
    treated <- c(treated1, treated2, treated3, treated4, treated5, treated6, treated7, treated8)
    
    # order data by all strata
    full_data <- full_data[order(full_data$all_strata), ]
    
    # Replicate treatment assignment for the number of people per cluster
    full_data$treatment <- rep(treated, each = n_units)
    
    linear_predictor <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment -
      (prop_effect[j])*full_data$proportion - 
      full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
    
    probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
    
    full_data$outcomes <- rbinom(n_clusters[j]*n_units, 1 , probability)
    
    full_data$within <- full_data$race - full_data$proportion
    
    full_data$within_cont <- full_data$continuous_pred - full_data$cont_ave
    
    
    # Now we define our models
    small_model <- update(small_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    strat_model <- update(strat_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    large_model <- update(large_model_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    kalbfleisch <- update(kalbfleisch_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    teerenstra <- update(teerenstra_1, newdata = full_data, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    
    
    
    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion + exp(continuous_pred) + cont_ave, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = binomial, corstr = "exchangeable"), error = function(e) NULL)
    
    singular_check[m,1] <- is.null(small_model)
    singular_check[m,2] <- is.null(strat_model)
    singular_check[m,3] <- is.null(large_model)
    singular_check[m,4] <- is.null(kalbfleisch)
    singular_check[m,5] <- is.null(teerenstra)
    
    if ((sum(singular_check[m,]) == 0) & (sum(is.null(gee_true) + is.null(gee_under)) == 0)){
      
      # small model
      marginal_model2 <- brmsmargins(
        small_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)], CI = 0.95)
      
      treatment_diff[m,2] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,3] <- sd(marginal_model2$Contrasts)
      
      # Repeat the same steps above for every model
      
      # true model
      marginal_model2 <- brmsmargins(
        large_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)], CI = 0.95)
      
      treatment_diff[m,4] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,5] <- sd(marginal_model2$Contrasts)
      
      # strat model
      marginal_model2 <- brmsmargins(
        strat_model,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)], CI = 0.95)
      
      treatment_diff[m,6] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,7] <- sd(marginal_model2$Contrasts)
      
      
      # kalbfleisch model
      marginal_model2 <- brmsmargins(
        kalbfleisch,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)], CI = 0.95)
      
      treatment_diff[m,8] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,9] <- sd(marginal_model2$Contrasts)
      
      
      # teerenstra model 
      marginal_model2 <- brmsmargins(
        teerenstra,
        at = data.frame(treatment = c(0, 1)),
        contrasts = matrix(c(-1, 1), nrow = 2),
        effects = "integrateoutRE", k = 100L, seed = seeds[((j-1)*3*n_sims + m + 2*n_sims)], CI = 0.95)
      
      treatment_diff[m,10] <- as.numeric(marginal_model2$ContrastSummary[1,1])
      treatment_diff[m,11] <- sd(marginal_model2$Contrasts)
      
      
      # gee true model
      marginal_model2 <- marginaleffects(gee_true)
      
      treatment_diff[m,12] <- summary(marginal_model2)[1,4]
      treatment_diff[m,13] <- summary(marginal_model2)[1,5]
      
      # gee under model
      marginal_model2 <- marginaleffects(gee_under)
      
      treatment_diff[m,14] <- summary(marginal_model2)[1,4]
      treatment_diff[m,15] <- summary(marginal_model2)[1,5]
      
    }
  }
  
  
  # predictors for the model
  cluster_mean <- rep(rgamma(18000, cluster_var[j], 1), each = 100)
  
  cluster <- rep(1:18000, each = 100)
  
  # extra_covariate <- rep(-1:(18000-2), each = 100)
  
  noise <- rnorm(18000*100, 0, 1)
  
  full_data <- data.frame(cluster_mean, cluster,   noise)
  
  names(full_data) <- c("cluster_mean", "cluster",   "noise")
  
  # proportion variable on which to stratify data
  prop <- runif(18000, 0, 1)
  
  # order the dataset by proportion variable
  full_data$proportion <- rep(prop, each = 100)
  full_data <- full_data[order(full_data$proportion), ]
  full_data$race <- rbinom(length(full_data$proportion), 1, full_data$proportion)
  
  # add strata identifier vector
  prop_strat <- c(1,2,3,4)
  prop_strata <- rep(prop_strat, each = ((18000/4)*100))
  
  ## Add strata vectors into the dataset
  full_data <- data.frame(full_data, prop_strata)
  
  ## Continuous predictor
  full_data$continuous_pred <- rnorm(18000*100, 0, 1)
  
  full_data$cont_ave <- ave(full_data$continuous_pred, full_data$cluster)
  
  # order the dataset by cluster mean variable
  full_data <- full_data[order(full_data$cont_ave), ]
  
  
  # # add continuous strata identifier vector
  cont_strat <- c(1,2)
  full_data$cont_strata <- rep(cont_strat, each = ((18000/2)*100))
  
  ## New Strata based on both cont and prop strata
  full_data$all_strata <- ifelse((full_data$prop_strata==1 & full_data$cont_strata==1), 1, 
                                 ifelse((full_data$prop_strata==1 & full_data$cont_strata==2), 2, 
                                        ifelse((full_data$prop_strata==2 & full_data$cont_strata==1), 3, 
                                               ifelse((full_data$prop_strata==2 & full_data$cont_strata==2), 4, 
                                                      ifelse((full_data$prop_strata==3 & full_data$cont_strata==1), 5, 
                                                             ifelse((full_data$prop_strata==3 & full_data$cont_strata==2), 6, 
                                                                    ifelse((full_data$prop_strata==4 & full_data$cont_strata==1), 7, 8))))) ))
  
  
  clust1 <- nrow(full_data[full_data$all_strata==1, ]) / 100
  clust2 <- nrow(full_data[full_data$all_strata==2, ]) / 100
  clust3 <- nrow(full_data[full_data$all_strata==3, ]) / 100
  clust4 <- nrow(full_data[full_data$all_strata==4, ]) / 100
  clust5 <- nrow(full_data[full_data$all_strata==5, ]) / 100
  clust6 <- nrow(full_data[full_data$all_strata==6, ]) / 100
  clust7 <- nrow(full_data[full_data$all_strata==7, ]) / 100
  clust8 <- nrow(full_data[full_data$all_strata==8, ]) / 100
  
  
  
  ## Stratify data by proportion race (here we use two strata)
  ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
  ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
  treated1 <- sample(c(rep(0, round(clust1/2)), rep(1, clust1 - round(clust1/2))), replace = F)
  treated2 <- sample(c(rep(0, round(clust2/2)), rep(1, clust2 - round(clust2/2))), replace = F)
  treated3 <- sample(c(rep(0, round(clust3/2)), rep(1, clust3 - round(clust3/2))), replace = F)
  treated4 <- sample(c(rep(0, round(clust4/2)), rep(1, clust4 - round(clust4/2))), replace = F)
  treated5 <- sample(c(rep(0, round(clust5/2)), rep(1, clust5 - round(clust5/2))), replace = F)
  treated6 <- sample(c(rep(0, round(clust6/2)), rep(1, clust6 - round(clust6/2))), replace = F)
  treated7 <- sample(c(rep(0, round(clust7/2)), rep(1, clust7 - round(clust7/2))), replace = F)
  treated8 <- sample(c(rep(0, round(clust8/2)), rep(1, clust8 - round(clust8/2))), replace = F)
  
  ## Combine the treatment vectors
  # full treatment vector
  treated <- c(treated1, treated2, treated3, treated4, treated5, treated6, treated7, treated8)
  
  # order data by all strata
  full_data <- full_data[order(full_data$all_strata), ]
  
  # Replicate treatment assignment for the number of people per cluster
  full_data$treatment <- rep(treated, each = 100)
  
  linear_predictor <- full_data$cluster_mean +
    (treat_effect[j])*full_data$treatment -
    (prop_effect[j])*full_data$proportion - 
    full_data$race + exp((full_data$continuous_pred / (cont_effect))) + full_data$cont_ave
  
  probability <- exp(linear_predictor) / (1 + exp(linear_predictor))
  
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  treatment_diff[,1] <- rep((mean(probability[which(full_data$treatment == 1)]) - mean(probability[which(full_data$treatment == 0)])), n_sims)
  
  # counter4later_total
  
  ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  treatment_diff <- na.omit(treatment_diff)
  
  ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  crt_design4_missp_RE_Gamma[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
                                              treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
                                              treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
                                              treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
                                              treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
                                              treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                              treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                              treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  crt_design4_missp_RE_Gamma[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_design4_missp_RE_Gamma[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  crt_design4_missp_RE_Gamma[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  
  crt_design4_missp_RE_Gamma[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  
  crt_design4_missp_RE_Gamma[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  
  crt_design4_missp_RE_Gamma[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  
  crt_design4_missp_RE_Gamma[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_design4_missp_RE_Gamma[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_design4_missp_RE_Gamma[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  crt_design4_missp_RE_Gamma[index, 22] <- mean(2*1.96*treatment_diff[,3])
  
  crt_design4_missp_RE_Gamma[index, 23] <- mean(2*1.96*treatment_diff[,7])
  
  crt_design4_missp_RE_Gamma[index, 24] <- mean(2*1.96*treatment_diff[,5])
  
  crt_design4_missp_RE_Gamma[index, 25] <- mean(2*1.96*treatment_diff[,9])
  
  crt_design4_missp_RE_Gamma[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_design4_missp_RE_Gamma[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_design4_missp_RE_Gamma[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  ############################################################################################################
  ############################################################################################################
  
  save(crt_design4_missp_RE_Gamma, file = "crt_design4_missp_RE_Gamma.Rdata")
  
  
  
}





