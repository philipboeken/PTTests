source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('helpers.R')
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))


# Input parameters
##############################################
n <- 300
m <- 500

err_sd <- 0.1
p_link <- 0.8

nonlin_options <- c(
  # linear,
  parabolic
  # sinusoidal
  # partial
  # circular
  # checkerboard
)


# Setup test
##############################################
get_data <- function(n, p_link, nonlin_options, err_sd) {
  link_nonlin <- rbinom(1, 1, p_link)
  
  res <- nonlin(nonlin_options, X=NULL, n=n)
  
  X <- res$X
  X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
  
  Y <- link_nonlin * res$Y
  Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
  
  return(list(X=X, Y=Y, label=as.numeric(!link_nonlin)))
  
}

get_results <- function(n, m){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n, p_link, nonlin_options, err_sd)
    return(data.frame(label=data$label,
                      cor=.pcor_wrapper(data$X, data$Y),
                      # splineGCM=.gcm_wrapper(data$X, data$Y),
                      # RCoT=.rcot_wrapper(data$X, data$Y),
                      # CCIT=.ccit_wrapper(data$X, data$Y),
                      Bayes=.bayes_wrapper(data$X, data$Y)))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- get_results(n, m)

stopCluster(.cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))
