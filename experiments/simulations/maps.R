# Maps
##############################################

### Interventions

mean_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (base + 2*theta)
}

variance_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (1+2*theta) * base
}

fixed_point <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (10^(-2*theta) * base + sd(base))
}

mixture <- function(base, C) {
  theta <- sample(1:3, 1)
  idx <- sample(c(-1,2), length(C), replace = TRUE)
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
  
  Y <- X_new
  
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
  
  Y <- X_new^2
  
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
  
  Y <- sin(6*2*pi*X_new / (max(X_new)-min(X_new)))
  
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
  
  b <- rbinom(length(X_new), 1, 0.2)
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

