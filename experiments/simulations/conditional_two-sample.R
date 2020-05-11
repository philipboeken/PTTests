source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('helpers.R')
library(foreach)
library(doParallel)


# Input parameters
##############################################
n <- 300
m <- 200

p_link <- 0.5

err_sd <- 0.1
nonlin_options <- c(
  linear,
  parabolic,
  sinusoidal,
  partial
)

p_two_sample <- 0.5
interv_options <- c(
  mean_shift,
  variance_shift,
  fixed_point,
  mixture
  # tails
)

p_ci <- 0.5


# Setup test
##############################################
get_data <- function(n, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options) {
  C <- rbinom(n, 1, p_two_sample)
  
  cond_indep <- rbinom(1, 1, p_ci)
  if (cond_indep) { # C -> Z -> X
    intervene <- rbinom(1, 1, p_link)
    Z <- intervene * do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
    
    link_nonlin <- rbinom(1, 1, p_link)
    X <- link_nonlin * nonlin(nonlin_options, Z)
    X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
  } else {
    if (runif(1) <= 0.5) { # C -> Z <- X
      X <- rnorm(n)
      
      link_nonlin <- rbinom(1, 1, p_link)
      Z <- link_nonlin * nonlin(nonlin_options, X)
      Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      
      intervene <- rbinom(1, 1, p_link)
      Z <- intervene * do_intervention(interv_options, Z, C) + (1-intervene) * Z
    } else { # C -> Z <- L -> X
      L <- rnorm(n)
      
      link_nonlin1 <- rbinom(1, 1, p_link)
      X <- link_nonlin1 * nonlin(nonlin_options, L)
      X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
      
      link_nonlin2 <- rbinom(1, 1, p_link)
      Z <- link_nonlin2 * L
      Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      
      intervene <- rbinom(1, 1, p_link)
      Z <- intervene * do_intervention(interv_options, Z, C) + (1-intervene) * Z
      
      link_nonlin <- link_nonlin1 & link_nonlin2
    }
  }
  
  cond_indep <- as.numeric(cond_indep | !link_nonlin | !intervene)
  
  return(list(C=C, Z=Z, X=X, label=cond_indep))
}

get_results <- function(n, m, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options)
    return(data.frame(
      label=data$label,
      pcor=.pcor_wrapper(data$C, data$X, data$Z),
      splineGCM=.gcm_wrapper(data$C, data$X, data$Z),
      RCoT=.rcot_wrapper(data$C, data$X, data$Z),
      CCIT=.ccit_wrapper(data$C, data$X, data$Z),
      Bayes=.bayes_wrapper(data$C, data$X, data$Z)
    ))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- get_results(n, m, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options)

stopCluster(.cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))
