source('indepTests/bayesUCITest.R')

bayes.CItest <- function (X, Y, Z = NULL, rho = 0.5, c = 1, max_depth = -1, qdist = qnorm, verbose = TRUE) {

  if (is.null(Z)) {
    if (verbose) {
      cat('Redirecting to UCI test\n')
    }
    return(bayes.UCItest(X, Y, c, max_depth, qdist, verbose))
  }

  max_depth <- ifelse(max_depth < 0, ceiling(log2(length(X))/2), max_depth)

  old_expressions <- options()$expressions
  options(expressions = max(2*max_depth, old_expressions))

  if (length(unique(Y)) == 2) {

    if (verbose) {
      cat('Performing two-sample test\n')
    }

    XYZ <- cbind(scale(X), Y, scale(Z))
    X1Z <- XYZ[XYZ[,2] == 0, c(1, 3)]
    X2Z <- XYZ[XYZ[,2] == 1, c(1, 3)]
    X <- X1Z[,1]
    Y <- X2Z[,1]

    p_x <- opt_pt(XYZ, 1, 3, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    p_x1 <- opt_pt(X1Z, 1, 2, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    p_x2 <- opt_pt(X2Z, 1, 2, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)

    bf <- p_x$p / (p_x1$p * p_x2$p) * 10^(p_x1$mag + p_x2$mag - p_x$mag)

  } else {

    X <- scale(X)
    Y <- scale(Y)
    XYZ <- cbind(X, Y, scale(Z))

    phi_x <- opt_pt(XYZ, 1, 3, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    phi_y <- opt_pt(XYZ, 2, 3, z_min=0, z_max=1, c=c, rho=rho, depth=1, max_depth, qdist)
    phi_xy <- opt_pt(XYZ, c(1, 2), 3, z_min=0, z_max=1, c=c, rho=rho, depth=1, max_depth, qdist)

    bf <- phi_x$p * phi_y$p / phi_xy$p * 10^(phi_xy$mag - phi_x$mag - phi_y$mag)

  }

  options(expressions = old_expressions)

  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf), X=X, Y=Y))
}

#########################################
# Helpers for OPT-PT

opt_pt <- function(data, target_idx, z_idx, z_min, z_max, c, rho, depth, max_depth, qdist) {
  local <- matrix(data[which(qdist(z_min) < data[, z_idx] & data[, z_idx] < qdist(z_max)),
                       target_idx], ncol=length(target_idx))
  logl <- pt_logl(local, low=rep(0, ncol(local)), up=rep(1, ncol(local)),
                  c=c, depth=1, max_depth, qdist)

  # if ((max_depth && depth == max_depth) || (!max_depth && nrow(local) <= 1)) {
  if (depth == max_depth) {
    mag <- floor(log10(exp(1))*logl)
    logl <- logl + log(10) * max(-mag - 5, 0)
    return(list(p=exp(logl), mag= max(-mag - 5, 0)))
  }


  phi_1 <- opt_pt(data, target_idx, z_idx, z_min, (z_min + z_max)/2,
                  c, rho, depth + 1, max_depth, qdist)
  phi_2 <- opt_pt(data, target_idx, z_idx, (z_min + z_max)/2, z_max,
                  c, rho, depth + 1, max_depth, qdist)

  return(apply_magnitude_correction(logl, phi_1$p, phi_1$mag, phi_2$p, phi_2$mag, rho))
}

apply_magnitude_correction <- function(p0, p1, mag1, p2, mag2, rho) {
  mag_p0 <- floor(log10(exp(1))*p0)
  mag <- mag1 + mag2
  if (exp(p0) != 0 && mag + mag_p0 > 300) {
    p <- rho * exp(p0)
    mag <- 0
  } else if (exp(p0) != 0){
    p <- rho * exp(log(10)*mag + p0) + (1 - rho) * p1 * p2
  } else {
    p <- (1 - rho) * p1 * p2
  }
  mag_p <- floor(log10(p))
  if (abs(mag_p) > 30) {
    p <- p * 10^(-mag_p - 5)
    mag <- mag + -mag_p - 5
  }
  return(list(p=p, mag=mag))
}

#########################################
# Helper for OPT-TPT

# opt_tpt <- function(XYZ, target_vars, z_min, z_max, c, rho, depth, max_depth) {
#   data <- matrix(XYZ[z_min < XYZ[,3] & XYZ[,3] < z_max, target_vars], ncol=length(target_vars))
#   likelihood <- tpt_likelihood(data, c, max_depth)
#
#   if (depth == max_depth || nrow(data) <= 1) {
#     return(likelihood)
#   }
#
#   phi_1 <- opt_tpt(XYZ, target_vars, z_min, (z_min + z_max)/2, c, rho, depth + 1, max_depth)
#   phi_2 <- opt_tpt(XYZ, target_vars, (z_min + z_max)/2, z_max, c, rho, depth + 1, max_depth)
#
#   return(rho * likelihood + (1-rho) * phi_1 * phi_2)
# }

