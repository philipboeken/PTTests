source('init.R', chdir=TRUE)

# Test wrappers
##############################################
.bayes_wrapper <- function (X, Y, Z=NULL) {
  return(bayes.CItest(X, Y, Z, verbose=FALSE)$p_H0)
}

.pcor_wrapper <- function(X, Y, Z=NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor.test(X, Y)$p.value)
  }
  return(pcor.test(X, Y, Z)$p.value)
}

.gcm_wrapper <- function(X, Y, Z=NULL) {
  if (is.null(Z)) {
    return(splineGCM(1, 2, c(), cbind(X, Y)))
  }
  return(splineGCM(1, 2, c(3), cbind(X, Y, Z)))
}

.rcot_wrapper <- function(X, Y, Z=NULL) {
  return(RCoT(X, Y, Z)$p)
}

.rcit_wrapper <- function(X, Y, Z=NULL) {
  return(RCIT(X, Y, Z)$p)
}

.ccit_wrapper <- function(X, Y, Z=NULL) {
  return(CCIT(X, Y, Z))
}


# Maps
##############################################

### Interventions

mean_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  theta <- 2
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

tails <- function(base, C) {
  # Fatter tails than base, when base is standard normal
  theta <- sample(1:3, 1)
  n <- length(C)
  (1-C) * base + C * rt(n, 10^theta)
}

log_mean_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  n <- length(C)
  (1-C) * rlnorm(n) + C * rlnorm(n, theta, 1)
}

log_variance_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  n <- length(C)
  (1-C) * rlnorm(n) + C * rlnorm(n, 0, 1+theta)
}

do_intervention <- function (int_options, base, C) {
  mapping <- sample(int_options, 1)[[1]]
  return(mapping(base, C))
}


### Nonlinearities

linear <- function(X=NULL, n=NULL) {
  X_new <- X
  if (is.null(X_new)) {
    X_new <- runif(n, 0, 2*pi)
  }
  
  Y <- 2 * X_new / 3
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}


parabolic <- function(X=NULL, n=NULL) {
  X_new <- X
  if (is.null(X_new)) {
    X_new <- runif(n, 0, 2*pi)
  }
  
  Y <- 2 * X_new^2 / 3
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}

sinusoidal <- function(X=NULL, n=NULL) {
  X_new <- X
  if (is.null(X_new)) {
    X_new <- runif(n, 0, 2*pi)
  }
  
  Y <- 2*sin(3*X_new)
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}

partial <- function(X=NULL, n=NULL) {
  X_new <- X
  if (is.null(X_new)) {
    X_new <- rnorm(n)
  }
  
  b <- rbinom(length(X_new), 1, 0.1)
  Y <- b*X_new + (1-b)*rnorm(length(X_new))
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}

circular <- function(X=NULL, n=NULL) {
  theta <- runif(n, 0, 2*pi)
  X_new <- cos(theta)
  Y <- sin(theta)
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}

checkerboard <- function(X=NULL, n=NULL) {
  i_x <- sample(c(1,2,3,4), n, replace=TRUE)
  i_y <- (i_x - rep(1, n)) %% 2 + sample(c(1,3), n, replace=TRUE)
  
  X_new <- (i_x + runif(n, 0, 1))
  Y <- (i_y + runif(n, 0, 1))
  
  if (is.null(X)) {
    return(list(X=X_new, Y=Y))
  }
  return(Y)
}

nonlin <- function(opts, X=NULL, n=NULL) {
  mapping <- sample(opts, 1)[[1]]
  return(mapping(X=X, n=n))
}

