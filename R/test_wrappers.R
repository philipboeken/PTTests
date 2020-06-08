# Test wrappers
##############################################
.polyatree_wrapper <- function (X, Y, Z = NULL) {
  polyatree_ci_test(X, Y, Z, verbose = FALSE)$p_H0
}

.ppcor_b_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(bayes_cor_test(X, Y)$p_H0)
  }
  bayes_pcor_test(X, Y, Z)$p_H0
}

.ppcor_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor.test(X, Y)$p.value)
  }
  ppcor::pcor.test(X, Y, Z)$p.value
}

.spcor_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor.test(X, Y, method = 'spearman')$p.value)
  }
  ppcor::pcor.test(X, Y, Z, method = 'spearman')$p.value
}

.gcm_wrapper <- function(X, Y, Z = NULL) {
  GeneralisedCovarianceMeasure::gcm.test(X, Y, Z, regr.method = "gam")$p.value
}

.rcot_wrapper <- function(X, Y, Z = NULL) {
  library(RCIT)
  RCoT(X, Y, Z)$p
}

.ccit_wrapper <- function(X, Y, Z = NULL) {
  .ccit <- reticulate::import('CCIT')
  .ccit <- .ccit$CCIT$CCIT
  if (is.null(Z) || length(Z) == 0) {
    return(.ccit(matrix(X, ncol = 1), matrix(Y, ncol = 1), NULL))
  }
  .ccit(matrix(X, ncol = 1), matrix(Y, ncol = 1), matrix(Z, ncol = 1))
}
