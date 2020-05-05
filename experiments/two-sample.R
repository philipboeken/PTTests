source('init.R', chdir=TRUE)
library(foreach)
library(doParallel)
library(distr)
library(ROCR)


# Define datasets
##############################################
gauss_mean_shift <- function(n, C, theta) {
  (1-C) * rnorm(n) + C * rnorm(n, theta, 1)
}

gauss_variance_shift <- function(n, C, theta) {
  (1-C) * rnorm(n) + C * rnorm(n, 0, 1+theta)
}

gauss_mixture <- function(n, C, theta) {
  (1-C) * rnorm(n) + C * r(UnivarMixingDistribution(Norm(mean=-theta, sd=1), 
                                                    Norm(mean=theta, sd=1),
                                                    mixCoeff=c(0.5, 0.5)))(n)
}

tails <- function(n, C, theta) {
  (1-C) * rnorm(n) + C * ((theta != 0) * rt(n, 10^theta) + (theta == 0) * rnorm(n))
}

log_mean_shift <- function(n, C, theta) {
  (1-C) * rlnorm(n) + C * rlnorm(n, theta, 1)
}

log_variance_shift <- function(n, C, theta) {
  (1-C) * rlnorm(n) + C * rlnorm(n, 0, 1+theta)
}


# Setup test
##############################################
get_data <- function(n) {
  theta <- sample(0:3, 1)
  p <- runif(1, 0.45, 0.65)
  C <- rbinom(n, 1, p)
  X <- sample(c(
    gauss_mean_shift,
    gauss_variance_shift,
    gauss_mixture,
    tails
    # log_mean_shift
    # log_variance_shift
  ), 1)[[1]](n, C, theta)
  
  data <- cbind(X, C)
  X1 <- data[data[,2] == 0, 1]
  X2 <- data[data[,2] == 1, 1]
  return(list(p=p, C=C, X=X, X1=X1, X2=X2, label=as.numeric(theta==0)))
  
}

get_results <- function(n, m){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n)
    return(data.frame(label=data$label,
             cor=cor.test(data$C, data$X)$p.value,
             # ks=ks.test(data$X1, data$X2)$p.value,
             # kruskal=kruskal.test(list(data$X1, data$X2))$p.value,
             # splineGCM=splineGCM(1, 2, c(), cbind(data$C, data$X)),
             RCoT=RCoT(data$C, data$X)$p,
             # RCIT=RCIT(data$C, data$X)$p,
             BayesTS=bayes.UCItest(data$X, data$C, verbose=FALSE)$p_H0))
  }
  return(result)
}


# Do test
##############################################
cores <- detectCores()
cl <- makeForkCluster(cores[1]-1)
registerDoParallel(cl)

results <- get_results(100, 200)

stopCluster(cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))


