source('independence_tests/splineGCM.R')
source('independence_tests/bayesCITest.R')
source('independence_tests/CCIT.R')
source('independence_tests/bayesPCor.R')

suppressWarnings(library(ppcor))
suppressWarnings(library(RCIT))
suppressWarnings(library(reticulate))
use_python('/usr/local/bin/python3')
.ccit <- import('CCIT')
.ccit <- .ccit$CCIT$CCIT

# Test wrappers
##############################################
.bayes_wrapper <- function (X, Y, Z=NULL) {
  return(bayes.CItest(X, Y, Z, verbose=FALSE)$p_H0)
}

.bcor_wrapper <- function(X, Y, Z=NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(bayes.cor.test(X, Y)$p_H0)
  }
  return(bayes.pcor.test(X, Y, Z)$p_H0)
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
