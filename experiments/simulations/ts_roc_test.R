source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('helpers.R')
library(foreach)
library(doParallel)


# Input parameters
##############################################
n <- 300
m <- 500

p_two_sample <- runif(1, 0.45, 0.65)
p_link <- 0.5

interv_options <- c(
  mean_shift
  # variance_shift
  # fixed_point,
  # mixture,
  # tails
  # log_mean_shift,
  # log_variance_shift
)


# Setup test
##############################################
get_data <- function(n, p_two_sample, p_link, interv_options) {
  C <- rbinom(n, 1, p_two_sample)
  intervene <- rbinom(1, 1, p_link)
  X <- intervene * do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
  
  return(list(C=C, X=X, label=as.numeric(!intervene)))
}

get_results <- function(n, m, p_two_sample, p_link, interv_options){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n, p_two_sample, p_link, interv_options)
    return(data.frame(
      label=data$label,
      cor=.pcor_wrapper(data$C, data$X),
      # splineGCM=.gcm_wrapper(data$C, data$X),
      # RCoT=.rcot_wrapper(data$C, data$X),
      # CCIT=.ccit_wrapper(data$C, data$X),
      Bayes=.bayes_wrapper(data$C, data$X)
    ))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- get_results(n, m, p_two_sample, p_link, interv_options)

stopCluster(.cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))

