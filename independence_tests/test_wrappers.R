source('independence_tests/polyatreeTests.R')
source('independence_tests/CCIT.R')
source('independence_tests/bayesPCor.R')

library(GeneralisedCovarianceMeasure)
library(ppcor)
library(RCIT)
library(reticulate)
use_python('/usr/local/bin/python3')
.ccit <- import('CCIT')
.ccit <- .ccit$CCIT$CCIT

# Test wrappers
##############################################
.polyatree_wrapper <- function (X, Y, Z = NULL) {
  return(polyatree_ci_test(X, Y, Z, verbose = FALSE)$p_H0)
}

.ppcor_b_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(bayes.cor.test(X, Y)$p_H0)
  }
  return(bayes.pcor.test(X, Y, Z)$p_H0)
}

.ppcor_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor.test(X, Y)$p.value)
  }
  return(pcor.test(X, Y, Z)$p.value)
}

.spcor_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor.test(X, Y, method = 'spearman')$p.value)
  }
  return(pcor.test(X, Y, Z, method = 'spearman')$p.value)
}

.gcm_wrapper <- function(X, Y, Z = NULL) {
  return(gcm.test(X, Y, Z, regr.method = "gam")$p.value)
}

.rcot_wrapper <- function(X, Y, Z = NULL) {
  return(RCoT(X, Y, Z)$p)
}

.rcit_wrapper <- function(X, Y, Z = NULL) {
  return(RCIT(X, Y, Z)$p)
}

.ccit_wrapper <- function(X, Y, Z = NULL) {
  return(CCIT(X, Y, Z))
}
