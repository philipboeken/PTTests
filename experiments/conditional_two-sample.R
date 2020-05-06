source('experiments/test_helpers.R')
library(foreach)
library(doParallel)


# Setup test
##############################################
n <- 300
m <- 200

p_link <- 1

err_sd <- 0.1
nonlin_options <- c(linear, parabolic, sinusoidal, partial)

p_two_sample <- runif(1, 0.45, 0.65)
interv_options <- c(mean_shift, variance_shift, mixture, tails)

p_ci <- 0.5


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
    if (runif(1) <= 0) { # C -> Z <- X
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
      Z <- link_nonlin2 * nonlin(nonlin_options, L)
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
    suffStat<-list(data=cbind(data$X, data$C, data$Z),
                   contextVars=c(2),verbose=FALSE,removeNAs=FALSE)
    return(data.frame(label=data$label,
                      pcor=pcor.test(data$X, data$C, data$Z)$p.value,
                      splineGCM=splineGCM(1, 2, c(3), cbind(data$X, data$C, data$Z)),
                      RCoT=RCoT(data$X, data$C, data$Z)$p,
                      # RCIT=RCIT(data$X, data$C, data$Z)$p,
                      bayes.CItest=bayes.CItest(data$X, data$C, data$Z, verbose=FALSE)$p_H0))
  }
  return(result)
}


# Do test
##############################################
cores <- detectCores()
cl <- makeForkCluster(cores[1]-1)
registerDoParallel(cl)

results <- get_results(n, m, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options)

stopCluster(cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))