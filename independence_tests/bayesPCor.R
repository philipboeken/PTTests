suppressMessages(library(invgamma))

bayes.cor.test <- function(X, Y) {
  bf <- .jzs_corbf(cor(X, Y), length(X))
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}

bayes.pcor.test <- function(X, Y, Z) {
  r0 <- sqrt(summary(lm(Y~Z))$r.squared)
  r1 <- sqrt(summary(lm(Y~Z+X))$r.squared)
  bf <- .jzs_partcorbf(r0, r1, 1, 2, length(X))
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}

# Main function to analytically calculate the BF for partial correlation
# see Wetzels, R. & Wagenmakers, E.-J. (2012). A default Bayesian 
# hypothesis test for correlations and partial correlations. Psychonomic 
# Bulletin & Review.
.jzs_partcorbf <- function(r0,r1,p0,p1,n) {
  int <- function(r,n,p,g) {
    a <- .5 * ((n-1-p) * log(1+g) - (n-1) * log(1+g * (1-r^2)))
    return(exp(a) * dinvgamma(g, shape=.5, scale=n/2))
  }
  
  bf10 <- integrate(int, lower=0, upper=25, r=r1, p=p1, n=n)$value /
    integrate(int, lower=0, upper=25, r=r0, p=p0, n=n)$value
  
  return(1/bf10)
}


# Main function to analytically calculate the BF for correlation
# see Wetzels, R. & Wagenmakers, E.-J. (2012). A default Bayesian 
# hypothesis test for correlations and partial correlations. Psychonomic 
# Bulletin & Review.
.jzs_corbf <- function(r,n) {
  int <- function(r,n,g) {
    a <- .5 * ((n-2) * log(1+g) - (n-1) * log(1+g * (1-r^2)))
    return(exp(a) * dinvgamma(g, shape=.5, scale=n/2))
  }
  
  bf10 <- integrate(int,lower=0,upper=25,r=r,n=n)$value
  
  return(1/bf10)
}