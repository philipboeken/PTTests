source('init.R', chdir=TRUE)
library(foreach)
library(doParallel)
library(ROCR)
library(mgcv)

# Define mappings
##############################################
mean_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (base + theta)
}

variance_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (1+theta) * base
}

mixture <- function(base, C) {
  theta <- sample(1:3, 1)
  idx <- sample(c(-1,1), length(C), replace = TRUE)
  (1-C) * base + C * (base + idx*theta)
}

do_intervention <- function (base, C) {
  mapping <- sample(c(
    mean_shift,
    variance_shift,
    mixture
  ), 1)[[1]]
  return(mapping(base, C))
}

linear <- function(X) {
  2*X / 3
}

parabolic <- function(X) {
  2*X^2 / 3
}

sinusoidal <- function(X) {
  2*sin(3*X)
}

partial <- function(X) {
  b <- rbinom(length(X), 1, 0.1)
  b*X + (1-b)*rnorm(length(X))
}

nonlin <- function(X) {
  mapping <- sample(c(
    linear,
    parabolic,
    sinusoidal,
    partial
  ), 1)[[1]]
  return(mapping(X))
}


# Setup test
##############################################
get_data <- function(n) {
  p <- runif(1, 0.45, 0.65)
  C <- rbinom(n, 1, p)
  
  cond_indep <- sample(c(0, 1), 1)
  if (cond_indep) { # C -> Z -> X
    intervene <- sample(c(1,1), 1) # NB: is always equal to 1
    Z <- intervene * do_intervention(rnorm(n), C) + (1-intervene) * rnorm(n)
    
    link_nonlin <- sample(c(1,1), 1) # NB: is always equal to 1
    errs <- rnorm(n, 0, runif(1, 1, 2.5))
    X <- link_nonlin * nonlin(Z) + errs
  } else {
    if (runif(1) <= 0) { # C -> Z <- X
      X <- rnorm(n)
      
      link_nonlin <- sample(c(1,1), 1) # NB: is always equal to 1
      errs <- rnorm(n, 0, runif(1, 1, 1.5))
      Z <- link_nonlin * nonlin(X) + errs
      
      intervene <- sample(c(1,1), 1) # NB: is always equal to 1
      Z <- intervene * do_intervention(Z, C) + (1-intervene) * Z
    } else { # C -> Z <- L -> X
      L <- rnorm(n)

      link_nonlin1 <- sample(c(1,1), 1) # NB: is always equal to 1
      errs <- rnorm(n, 0, runif(1, 1, 2.5))
      X <- link_nonlin1 * nonlin(L) + errs

      link_nonlin2 <- sample(c(1,1), 1) # NB: is always equal to 1
      errs <- rnorm(n, 0, runif(1, 1, 2.5))
      Z <- link_nonlin2 * nonlin(L) + errs

      intervene <- sample(c(1,1), 1) # NB: is always equal to 1
      Z <- intervene * do_intervention(Z, C) + (1-intervene) * Z

      link_nonlin <- link_nonlin1 & link_nonlin2
    }
  }
  
  cond_indep <- as.numeric(cond_indep | !link_nonlin | !intervene)
  
  return(list(C=C, Z=Z, X=X, label=cond_indep))
}

get_results <- function(n, m){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n)
    suffStat<-list(data=cbind(data$X, data$C, data$Z),
                   contextVars=c(2),verbose=FALSE,removeNAs=FALSE)
    return(data.frame(label=data$label,
                      # pcor=pcor.test(data$X, data$C, data$Z)$p.value,
                      # gaussCIsincontest=gaussCIsincontest(1, 2, c(3), suffStat),
                      # gaussCIcontexttest=gaussCIcontexttest(1, 2, c(3), suffStat),
                      # splineGCM=splineGCM(1, 2, c(3), cbind(data$X, data$C, data$Z)),
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

results <- get_results(600, 400)

stopCluster(cl)


# Process results
##############################################
plot(pplot_roc(results[,1], results[,-1]))