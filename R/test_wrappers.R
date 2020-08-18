# Test wrappers
##############################################
.pt_wrapper <- function(X, Y, Z = NULL) {
  pt_ci_test(X, Y, Z, verbose = FALSE)$p_H0
}

.pt_wrapper_continuous <- function(X, Y, Z = NULL) {
  if (is.null(Z)) {
    return(pt_independence_test(X, Y, verbose = FALSE)$p_H0)
  }

  return(pt_continuous_ci_test(X, Y, Z)$p_H0)
}

.ppcor_b_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(cor_test_bayesian(X, Y)$p_H0)
  }
  pcor_test_bayesian(X, Y, Z)$p_H0
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
  if (!"RCIT" %in% (.packages()))
    library(RCIT)
  RCoT(X, Y, Z)$p
}

.kcit_wrapper <- function(X, Y, Z = NULL) {
  if (!"RCIT" %in% (.packages()))
    library(RCIT)
  KCIT(X, Y, Z)
}

.ccit_wrapper <- function(X, Y, Z = NULL) {
  .ccit <- reticulate::import('CCIT')
  .ccit <- .ccit$CCIT$CCIT
  if (is.null(Z) || length(Z) == 0) {
    return(reticulate::py_suppress_warnings(.ccit(matrix(X, ncol = 1), matrix(Y, ncol = 1), NULL)))
  }
  reticulate::py_suppress_warnings(.ccit(matrix(X, ncol = 1), matrix(Y, ncol = 1), matrix(Z, ncol = 1)))
}

.kruskal_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(kruskal.test(Y, X)$p.value)
  }
  .ppcor_wrapper(X, Y, Z)
}

.kolmogorov_wrapper <- function(X, Y, Z = NULL) {
  if ((is.null(Z) || length(Z) == 0) && .is_discrete(X)) {
    data <- cbind(X, Y)
    X1 <- data[data[, 1] == 0, 2]
    X2 <- data[data[, 1] == 1, 2]
    return(ks.test(X1, X2)$p.value)
  }
  .ppcor_wrapper(X, Y, Z)
}

.wilcoxon_wrapper <- function(X, Y, Z = NULL) {
  if (is.null(Z) || length(Z) == 0) {
    return(wilcox.test(Y, X)$p.value)
  }
  .ppcor_wrapper(X, Y, Z)
}
