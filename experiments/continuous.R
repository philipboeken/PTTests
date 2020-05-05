source('init.R', chdir=TRUE)
library(foreach)
library(doParallel)
library(ROCR)

# Define datasets
##############################################
linear <- function(n, sigma, indep) {
  X <- rnorm(n)
  Y <- (1-indep)*2*X / 3 + rnorm(n, 0, sigma)
  list(X=X, Y=Y)
}

parabolic <- function(n, sigma, indep) {
  X <- rnorm(n)
  Y <- (1-indep)*2*X^2 / 3 + rnorm(n, 0, sigma)
  list(X=X, Y=Y)
}

sinusoidal <- function(n, sigma, indep) {
  X <- runif(n, 0, 2*pi)
  Y <- (1-indep)*2*sin(3*X) + rnorm(n, 0, sigma)
  list(X=X, Y=Y)
}

circular <- function(n, sigma, indep) {
  theta <- runif(n, 0, 2*pi)
  X <- (1-indep)*10*cos(theta) + rnorm(n, 0, sigma)
  Y <- 10*sin(theta) + rnorm(n, 0, sigma)
  list(X=X, Y=Y)
}

checkerboard <- function(n, sigma, indep) {
  i_x <- sample(c(1,2,3,4), n, replace=TRUE)
  i_y <- (i_x - rep(1, n)) %% 2 + sample(c(1,3), n, replace=TRUE)
  X <- (1-indep)*10*(i_x + runif(n, 0, 1)) + rnorm(n, 0, sigma)
  Y <- 10*(i_y + runif(n, 0, 1)) + rnorm(n, 0, sigma)
  list(X=X, Y=Y)
}


# Setup test
##############################################
get_data <- function(n) {
  sigma <- runif(1, 1, 2.5)
  indep <- sample(c(0, 1), 1)
  res <- sample(c(
    linear,
    parabolic,
    sinusoidal,
    circular,
    checkerboard
  ), 1)[[1]](n, sigma, indep)
  
  return(list(X=res$X, Y=res$Y, true=indep))
  
}

get_results <- function(n, m){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n)
    return(data.frame(true=data$true,
                      cor=cor.test(data$X, data$Y)$p.value,
                      splineGCM=splineGCM(1, 2, c(), cbind(data$X, data$Y)),
                      RCoT=RCoT(data$X, data$Y)$p,
                      # RCIT=RCIT(data$X, data$Y)$p,
                      Bayes=bayes.UCItest(data$X, data$Y, verbose=FALSE)$p_H0))
  }
  return(result)
}


# Do test
##############################################
cores <- detectCores()
cl <- makeForkCluster(cores[1]-1)
registerDoParallel(cl)

results <- get_results(200, 200)

stopCluster(cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))
