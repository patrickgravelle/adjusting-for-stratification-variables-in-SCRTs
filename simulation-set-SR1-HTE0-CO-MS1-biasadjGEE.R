
library(dplyr)
library(geepack)
library(rstan)
library(brms)
library(brmsmargins)
library(marginaleffects)

library(geeCRT)
library(geesmv)

library(gee)
library(fBasics)

gee_var_md <- function (formula, id, family = gaussian, data, corstr = "independence") 
{
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  gee.fit <- gee(formula, data = data, id = id, family = family, 
                 corstr = corstr)
  beta_est <- gee.fit$coefficient
  alpha <- gee.fit$working.correlation[1, 2]
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- gee.fit$id
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster), 
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster, 
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation, 
    )
  }
  else {
    print(corstr)
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian", binomial = "binomial", 
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      print(family)
      stop("'family' not recognized")
    }
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)
  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])], 
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "gaussian") {
      xx <- t(covariate) %*% solve(var_i) %*% covariate
      step11 <- step11 + xx
    }
    else if (family == "poisson") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))), 
                 cluster$n[i]) %*% var_i %*% diag(sqrt(c(exp(covariate %*% 
                                                               beta_est))), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
    else if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 + 
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + 
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*% 
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 + 
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    }
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])], 
                                  mat$subj == unique(data$id)[i]))
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if (family == "gaussian") {
      xy <- t(covariate) %*% solve(var_i) %*% solve(cormax.ind(cluster$n[i]) - 
                                                      covariate %*% solve(step11) %*% t(covariate) %*% 
                                                      solve(var_i)) %*% (y - covariate %*% beta_est)
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
    else if (family == "poisson") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))), 
                 cluster$n[i]) %*% var_i %*% diag(sqrt(c(exp(covariate %*% 
                                                               beta_est))), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% solve(cormax.ind(cluster$n[i]) - 
                                           D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*% 
        (y - exp(covariate %*% beta_est))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
    else if (family == "binomial") {
      D <- mat.prod(covariate, exp(covariate %*% beta_est)/((1 + 
                                                               exp(covariate %*% beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + 
                                                       exp(covariate %*% beta_est))^2)), cluster$n[i]) %*% 
        var_i %*% diag(sqrt(c(exp(covariate %*% beta_est)/(1 + 
                                                             exp(covariate %*% beta_est))^2)), cluster$n[i])
      xy <- t(D) %*% solve(Vi) %*% solve(cormax.ind(cluster$n[i]) - 
                                           D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*% 
        (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*% 
                                                    beta_est)))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + vec(xy %*% t(xy))
      p[, i] <- vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*% 
    step14 %*% kronecker(solve(step11), solve(step11))
  return(list(cov.beta = cov.beta, cov.var = cov.var))
}


## Coverage, MSE, and Bias matrix
crt_updated_APR29_con_gee <- matrix(NA, nrow = 243, ncol = 28)

# counter4later_matrix

## Indices for placing values properly
index <- 0


# what's varying
# strata (but this is dealt with in the code)
# number of clusters
# treatment effect
# prop effect
# cluster variance

# these will need to be the full length vectors rather just unique values

# crt_updated_APR29_con_gee should be fine because index <- index + 1 should be fine through the loop


n_units <- 20
n_sims <- 100


# 2 strata 


# ## Define the number of clusters, units, and simulations
# n_clusters <- rep(c(24,24,24,36,36,36,48,48,48), 27)
# 
# ## Define the effect of the treatment
# treat_effect <- rep(rep(c(0.25,1,4), each = 27), 3)
# 
# ## Define the effect of the proportion of race
# prop_effect <- rep(c(1,2,3), each = 81)
# 
# cluster_var <- rep(c(rep(1/3, 9), rep(1, 9), rep(3, 9)), 9)




## Define the number of clusters, units, and simulations
n_clusters <- rep(c(24,36,48), 27)

## Define the effect of the treatment
treat_effect <- rep(rep(c(0.25,1,4), each = 9), 3)

## Define the effect of the proportion of race
prop_effect <- rep(c(1,2,3), each = 27)

cluster_var <- rep(c(rep(sqrt(1/3), 3), rep(1, 3), rep(sqrt(3), 3)), 9)


seeds <- (1:(length(cluster_var)*n_sims*3) + (length(cluster_var)*n_sims*3)*47)

priorset <- c(set_prior("normal(0,5)", class = "b"),
              set_prior("student_t(3,0,5)", class = "sd"))


# Broke at iteration 46 so let's try again from here
# for (j in 1:length(cluster_var)){
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
    cluster_mean <- rep(rnorm(n_clusters[j], 0, cluster_var[j]), each = n_units)
    
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
    
    ## Stratify data by proportion race (here we use two strata)
    ## So take the bottom 50% of clusters in terms of lower proportion and randomly allocate treatment to half of this strata
    ### This is just a vector of 0s and 1s that we expand by the number of people per cluster later
    treated1 <- sample(c(rep(0, round(n_clusters[j]/4)), rep(1, (n_clusters[j]/2 -round(n_clusters[j]/4)))), replace = F)
    
    ## Another treatment vector for the upper 50% of cluster proportions
    treated2 <- sample(c(rep(0, round(n_clusters[j]/4)), rep(1, (n_clusters[j]/2 -round(n_clusters[j]/4)))), replace = F)
    
    ## Combine the treatment vectors
    # full treatment vector
    treated <- c(treated1, treated2)
    # Replicate treatment assignment for the number of people per cluster
    treatment <- rep(treated, each = n_units)
    
    # add strata identifier vector
    strat <- c(1,2)
    strata <- rep(strat, each = ((n_clusters[j]/2)*n_units))
    
    ## Add strata and treatment vectors into the dataset
    full_data <- data.frame(full_data, treatment, strata)
    
    ## Generate our outcome by taking the sum of the following:
    # Cluster mean
    # treatment assignment (0 or 1) times the defined treatment effect
    # proportion race (cluster level) times the defined proportion effect
    # race (0 or 1)
    # noise 
    full_data$outcomes <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment +
      (prop_effect[j])*full_data$proportion + 
      full_data$noise + full_data$race
    
    ## For one of the literature models we create the "within cluster" correlation variable defined as follows
    ## individual value of race divided by the mean (here I used the proportion rather than the specific cluster mean)
    ## This would only effect this model's specific results
    full_data$within <- full_data$race - full_data$proportion
    
    # Now we define our models

    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion, 
                                id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), 
                          error = function(e) NULL)
    
    
    gee_true_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race + proportion, 
                                       id = "cluster", data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    gee_under_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race, id = "cluster", 
                                        data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    
    
    if ((sum(is.null(gee_true) + is.null(gee_under) + is.null(gee_true_se) + is.null(gee_under_se)) == 0)){
      
      
      # gee true model
      treatment_diff[m,12] <- summary(gee_true)$coefficients[2,1]
      treatment_diff[m,13] <- sqrt(gee_true_se$cov.beta[2,2])
      
      # gee under model
      treatment_diff[m,14] <- summary(gee_under)$coefficients[2,1]
      treatment_diff[m,15] <- sqrt(gee_under_se$cov.beta[2,2])
      
    }
    
    
  }
  
  
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  # treatment_diff[,1] <- rep(treat_effect[j], n_sims)
  
  # counter4later_total
  
  trt_diff_test <- na.omit(treatment_diff[,c(12:15)])
  
  treatment_diff <- matrix(0, nrow = nrow(trt_diff_test), ncol = 15)
  
  treatment_diff[,1] <- rep(treat_effect[j], nrow(trt_diff_test))
  
  treatment_diff[,c(12:15)] <- trt_diff_test[,c(1:4)]
  
  
  # ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  # treatment_diff <- na.omit(treatment_diff)
  # 
  # ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  # crt_updated_APR29_con_gee[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                             treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                             treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  # crt_updated_APR29_con_gee[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  # crt_updated_APR29_con_gee[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  # 
  # crt_updated_APR29_con_gee[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  # 
  # crt_updated_APR29_con_gee[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  # 
  # crt_updated_APR29_con_gee[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  # 
  # crt_updated_APR29_con_gee[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_updated_APR29_con_gee[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_updated_APR29_con_gee[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  # crt_updated_APR29_con_gee[index, 22] <- mean(2*1.96*treatment_diff[,3])
  # 
  # crt_updated_APR29_con_gee[index, 23] <- mean(2*1.96*treatment_diff[,7])
  # 
  # crt_updated_APR29_con_gee[index, 24] <- mean(2*1.96*treatment_diff[,5])
  # 
  # crt_updated_APR29_con_gee[index, 25] <- mean(2*1.96*treatment_diff[,9])
  # 
  # crt_updated_APR29_con_gee[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_updated_APR29_con_gee[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_updated_APR29_con_gee[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  
  
  
  
  
  
  
  
  
  
  
  # 3 Strata
  
  
  treatment_diff <- matrix(NA, nrow = n_sims, ncol = 15)
  
  index <- index + 1
  
  singular_check <- matrix(NA, nrow = n_sims, ncol = 5)
  
  
  ################################
  
  for (m in 1:n_sims){
    
    set.seed(seeds[((j-1)*3*n_sims + m + 1*n_sims)])
    
    
    # predictors for the model
    cluster_mean <- rep(rnorm(n_clusters[j], 0, cluster_var[j]), each = n_units)
    
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
    
    # treatment assignment for each strata
    treated1 <- sample(c(rep(0, round(n_clusters[j]/6)), rep(1, (n_clusters[j]/3 -round(n_clusters[j]/6)))), replace = F)
    
    treated2 <- sample(c(rep(0, round(n_clusters[j]/6)), rep(1, (n_clusters[j]/3 -round(n_clusters[j]/6)))), replace = F)
    
    treated3 <- sample(c(rep(0, round(n_clusters[j]/6)), rep(1, (n_clusters[j]/3 -round(n_clusters[j]/6)))), replace = F)
    
    
    # full treatment vector
    treated <- c(treated1, treated2, treated3)
    treatment <- rep(treated, each = n_units)
    
    # add strata identifier
    strat <- c(1,2,3)
    strata <- rep(strat, each = ((n_clusters[j]/3)*n_units))
    
    full_data <- data.frame(full_data, treatment, strata)
    
    full_data$outcomes <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment +
      (prop_effect[j])*full_data$proportion +
      full_data$noise + full_data$race
    
    full_data$within <- full_data$race - full_data$proportion
    
    
    
    
    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion, 
                                id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), 
                          error = function(e) NULL)
    
    
    gee_true_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race + proportion, 
                                       id = "cluster", data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    gee_under_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race, id = "cluster", 
                                        data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    
    
    if ((sum(is.null(gee_true) + is.null(gee_under) + is.null(gee_true_se) + is.null(gee_under_se)) == 0)){
      
      
      # gee true model
      treatment_diff[m,12] <- summary(gee_true)$coefficients[2,1]
      treatment_diff[m,13] <- sqrt(gee_true_se$cov.beta[2,2])
      
      # gee under model
      treatment_diff[m,14] <- summary(gee_under)$coefficients[2,1]
      treatment_diff[m,15] <- sqrt(gee_under_se$cov.beta[2,2])
      
    }
    
  }
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  # treatment_diff[,1] <- rep(treat_effect[j], n_sims)
  
  trt_diff_test <- na.omit(treatment_diff[,c(12:15)])
  
  treatment_diff <- matrix(0, nrow = nrow(trt_diff_test), ncol = 15)
  
  treatment_diff[,1] <- rep(treat_effect[j], nrow(trt_diff_test))
  
  treatment_diff[,c(12:15)] <- trt_diff_test[,c(1:4)]
  
  # counter4later_total
  
  # ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  # treatment_diff <- na.omit(treatment_diff)
  # 
  # ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  # crt_updated_APR29_con_gee[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                                      treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                                      treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  # crt_updated_APR29_con_gee[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  # crt_updated_APR29_con_gee[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  # 
  # crt_updated_APR29_con_gee[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  # 
  # crt_updated_APR29_con_gee[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  # 
  # crt_updated_APR29_con_gee[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  # 
  # crt_updated_APR29_con_gee[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_updated_APR29_con_gee[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_updated_APR29_con_gee[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  # crt_updated_APR29_con_gee[index, 22] <- mean(2*1.96*treatment_diff[,3])
  # 
  # crt_updated_APR29_con_gee[index, 23] <- mean(2*1.96*treatment_diff[,7])
  # 
  # crt_updated_APR29_con_gee[index, 24] <- mean(2*1.96*treatment_diff[,5])
  # 
  # crt_updated_APR29_con_gee[index, 25] <- mean(2*1.96*treatment_diff[,9])
  # 
  # crt_updated_APR29_con_gee[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_updated_APR29_con_gee[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_updated_APR29_con_gee[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  
  
  # 4 Strata
  
  
  treatment_diff <- matrix(NA, nrow = n_sims, ncol = 15)
  
  index <- index + 1
  
  singular_check <- matrix(NA, nrow = n_sims, ncol = 5)
  
  
  ################################
  
  for (m in 1:n_sims){
    
    set.seed(seeds[((j-1)*3*n_sims + m + 2*n_sims)])
    
    
    # predictors for the model
    cluster_mean <- rep(rnorm(n_clusters[j], 0, cluster_var[j]), each = n_units)
    
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
    # treatment assignment for each strata
    treated1 <- sample(c(rep(0, round(n_clusters[j]/8)), rep(1, (n_clusters[j]/4 -round(n_clusters[j]/8)))), replace = F)
    
    treated2 <- sample(c(rep(0, (n_clusters[j]/4 -round(n_clusters[j]/8))), rep(1, round(n_clusters[j]/8))), replace = F)
    
    treated3 <- sample(c(rep(0, round(n_clusters[j]/8)), rep(1, (n_clusters[j]/4 -round(n_clusters[j]/8)))), replace = F)
    
    treated4 <- sample(c(rep(0, (n_clusters[j]/4 -round(n_clusters[j]/8))), rep(1, round(n_clusters[j]/8))), replace = F)
    
    
    # full treatment vector
    treated <- c(treated1, treated2, treated3, treated4)
    treatment <- rep(treated, each = n_units)
    
    # add strata identifier
    strat <- c(1,2,3,4)
    strata <- rep(strat, each = ((n_clusters[j]/4)*n_units))
    
    full_data <- data.frame(full_data, treatment, strata)
    
    full_data$outcomes <- full_data$cluster_mean +
      (treat_effect[j])*full_data$treatment +
      (prop_effect[j])*full_data$proportion +
      full_data$noise + full_data$race
    
    full_data$within <- full_data$race - full_data$proportion
    
    
    
    
    gee_true <- tryCatch(geeglm(outcomes ~ factor(treatment) + race + proportion, 
                                id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    gee_under <- tryCatch(geeglm(outcomes ~ factor(treatment) + race, id = cluster, data = full_data, family = gaussian, corstr = "exchangeable"), 
                          error = function(e) NULL)
    
    
    gee_true_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race + proportion, 
                                       id = "cluster", data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    gee_under_se <- tryCatch(gee_var_md(outcomes ~ factor(treatment) + race, id = "cluster", 
                                        data = full_data, family = gaussian, corstr = "exchangeable"), error = function(e) NULL)
    
    
    
    
    if ((sum(is.null(gee_true) + is.null(gee_under) + is.null(gee_true_se) + is.null(gee_under_se)) == 0)){
      
      
      # gee true model
      treatment_diff[m,12] <- summary(gee_true)$coefficients[2,1]
      treatment_diff[m,13] <- sqrt(gee_true_se$cov.beta[2,2])
      
      # gee under model
      treatment_diff[m,14] <- summary(gee_under)$coefficients[2,1]
      treatment_diff[m,15] <- sqrt(gee_under_se$cov.beta[2,2])
      
    }
    
  }
  
  ## Add a vector of the true treatment effect to make confidence interval computation easier
  # treatment_diff[,1] <- rep(treat_effect[j], n_sims)
  
  trt_diff_test <- na.omit(treatment_diff[,c(12:15)])
  
  treatment_diff <- matrix(0, nrow = nrow(trt_diff_test), ncol = 15)
  
  treatment_diff[,1] <- rep(treat_effect[j], nrow(trt_diff_test))
  
  treatment_diff[,c(12:15)] <- trt_diff_test[,c(1:4)]
  
  # counter4later_total
  
  # ## remove rows with any NAs (the only NAs that appear would be from us skipping an iteration due to model(s) with a singular fit)
  # treatment_diff <- na.omit(treatment_diff)
  # 
  # ## using ifelse function to compute coverage by summing over TRUE's for if the true treatment effect is within the 95% Confidence Interval
  # crt_updated_APR29_con_gee[index, 1] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,2] - 1.96*treatment_diff[,3]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,2] + 1.96*treatment_diff[,3]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 2] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,6] - 1.96*treatment_diff[,7]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,6] + 1.96*treatment_diff[,7]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 3] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,4] - 1.96*treatment_diff[,5]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,4] + 1.96*treatment_diff[,5]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 4] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,8] - 1.96*treatment_diff[,9]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,8] + 1.96*treatment_diff[,9]), 1, 0)) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 5] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,10] - 1.96*treatment_diff[,11]) & 
  #                                            treatment_diff[,1] <= (treatment_diff[,10] + 1.96*treatment_diff[,11]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 6] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,12] - 1.96*treatment_diff[,13]) & 
                                                      treatment_diff[,1] <= (treatment_diff[,12] + 1.96*treatment_diff[,13]), 1, 0)) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 7] <- sum(ifelse(treatment_diff[,1] >= (treatment_diff[,14] - 1.96*treatment_diff[,15]) & 
                                                      treatment_diff[,1] <= (treatment_diff[,14] + 1.96*treatment_diff[,15]), 1, 0)) / nrow(treatment_diff)
  
  # MSE of every model
  # crt_updated_APR29_con_gee[index, 8] <- (sum(((treatment_diff[,2] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 9] <- (sum(((treatment_diff[,6] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 10] <- (sum(((treatment_diff[,4] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 11] <- (sum(((treatment_diff[,8] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  # 
  # crt_updated_APR29_con_gee[index, 12] <- (sum(((treatment_diff[,10] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 13] <- (sum(((treatment_diff[,12] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  crt_updated_APR29_con_gee[index, 14] <- (sum(((treatment_diff[,14] - treatment_diff[,1])^2))) / nrow(treatment_diff)
  
  # summed absolute bias of every model
  # crt_updated_APR29_con_gee[index, 15] <- mean(abs(treatment_diff[,1] - treatment_diff[,2]))
  # 
  # crt_updated_APR29_con_gee[index, 16] <- mean(abs(treatment_diff[,1] - treatment_diff[,6]))
  # 
  # crt_updated_APR29_con_gee[index, 17] <- mean(abs(treatment_diff[,1] - treatment_diff[,4]))
  # 
  # crt_updated_APR29_con_gee[index, 18] <- mean(abs(treatment_diff[,1] - treatment_diff[,8]))
  # 
  # crt_updated_APR29_con_gee[index, 19] <- mean(abs(treatment_diff[,1] - treatment_diff[,10]))
  
  crt_updated_APR29_con_gee[index, 20] <- mean(abs(treatment_diff[,1] - treatment_diff[,12]))
  
  crt_updated_APR29_con_gee[index, 21] <- mean(abs(treatment_diff[,1] - treatment_diff[,14]))
  
  # Average Interval width
  
  # crt_updated_APR29_con_gee[index, 22] <- mean(2*1.96*treatment_diff[,3])
  # 
  # crt_updated_APR29_con_gee[index, 23] <- mean(2*1.96*treatment_diff[,7])
  # 
  # crt_updated_APR29_con_gee[index, 24] <- mean(2*1.96*treatment_diff[,5])
  # 
  # crt_updated_APR29_con_gee[index, 25] <- mean(2*1.96*treatment_diff[,9])
  # 
  # crt_updated_APR29_con_gee[index, 26] <- mean(2*1.96*treatment_diff[,11])
  
  crt_updated_APR29_con_gee[index, 27] <- mean(2*1.96*treatment_diff[,13])
  
  crt_updated_APR29_con_gee[index, 28] <- mean(2*1.96*treatment_diff[,15])
  
  
  ############################################################################################################
  ############################################################################################################
  
  save(crt_updated_APR29_con_gee, file = "crt_updated_APR29_con_gee.Rdata")
  
}







