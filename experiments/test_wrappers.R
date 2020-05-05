source('init.R', chdir=TRUE)

bayes_ci_wrapper <- function (X, Y, Z=NULL) {
  bayes.CItest(X, Y, Z, verbose=FALSE)$p_H0
}

pcor_wrapper <- function(X, Y, Z=NULL) {
  if (is.null(Z)) {
    cor.test(X, Y)$p.value
  }
  pcor.test(X, Y, Z)$p.value
}

gcm_wrapper <- function(X, Y, Z=NULL) {
  if (is.null(Z)) {
    splineGCM(1, 2, c(), cbind(X, Y))
  }
  splineGCM(1, 2, c(3), cbind(X, Y, Z))
}

rcot_wrapper <- function(X, Y, Z=NULL) {
  RCoT(X, Y, Z)$p
}

rcit_wrapper <- function(X, Y, Z=NULL) {
  RCIT(X, Y, Z)$p
}
