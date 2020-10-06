pt_ci_test <- function(X, Y, Z = NULL, rho = 0.5, c = 1,
                       max_depth = -1, qdist = qnorm,
                       verbose = TRUE, log_BF = FALSE) {
  if (is.null(Z)) {
    if (verbose) cat('Redirecting to independence test\n')
    return(pt_independence_test(X, Y, c, max_depth, qdist, verbose, log_BF))
  }

  if (.is_discrete(X) || .is_discrete(Y)) {
    if (verbose) cat('Performing conditional two-sample test\n')
    return(pt_d_sample_ci_test(X, Y, Z, rho, c, max_depth, qdist, log_BF))
  }

  if (verbose) cat('Performing continuous conditional independence test\n')
  return(pt_continuous_ci_test(X, Y, Z, rho, c, max_depth, qdist, log_BF))
}

pt_d_sample_ci_test <- function(X, Y, Z, rho = 0.5, c = 1,
                                max_depth = -1, qdist = qnorm, log_BF = FALSE) {
  old_expressions <- options()$expressions
  options(expressions = max(max_depth, old_expressions))

  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X)) / 2)), max_depth)

  if (.is_discrete(X) && .is_discrete(Y) || length(X) <= 2) {
    return(list(bf = 1, p_H0 = 1 / 2, p_H1 = 1 / 2))
  }

  binary <- if (.is_discrete(X)) X else Y
  continuous <- if (.is_discrete(X)) Y else X

  Z <- matrix(Z, nrow = length(X))
  XYZ <- cbind(scale(continuous), binary, scale(Z))

  p_H0 <- .condopt_marginal_likelihood(XYZ, target_idx = 1, z_idx = 3,
                                       z_min = 0, z_max = 1,
                                       c = c, rho = rho, depth = 1, max_depth, qdist)

  discrete_values <- if (length(unique(binary)) == 2) binary[1] else unique(binary)

  p_H1 <- max(sapply(discrete_values, function(i) {
    X1Z <- matrix(XYZ[XYZ[, 2] == i,], ncol = ncol(XYZ))
    X2Z <- matrix(XYZ[XYZ[, 2] != i,], ncol = ncol(XYZ))

    p_x1 <- .condopt_marginal_likelihood(X1Z, target_idx = 1, z_idx = 3,
                                         z_min = 0, z_max = 1,
                                         c = c, rho = rho, depth = 1, max_depth, qdist)
    p_x2 <- .condopt_marginal_likelihood(X2Z, target_idx = 1, z_idx = 3,
                                         z_min = 0, z_max = 1,
                                         c = c, rho = rho, depth = 1, max_depth, qdist)
    p_x1 + p_x2
  }))

  n_hypotheses <- length(discrete_values)
  bf <- p_H0 - p_H1 + log(n_hypotheses) # Bayes Factor with a Bonferroni-type correction

  options(expressions = old_expressions)

  if (log_BF)
    return(list(bf = bf, p_H0 = NULL, p_H1 = NULL))

  bf <- exp(bf)
  return(list(bf = bf, p_H0 = 1 - 1 / (1 + bf), p_H1 = 1 / (1 + bf)))
}

# Implements a continuous conditional independence test.
# See Teymur, O. & Filippi, S. (2020). A Bayesian 
# nonparametric test for conditional independence.
# https://arxiv.org/abs/1910.11219
pt_continuous_ci_test <- function(X, Y, Z = NULL, rho = 0.5, c = 1,
                                  max_depth = -1, qdist = qnorm, log_BF = FALSE) {
  old_expressions <- options()$expressions
  options(expressions = max(max_depth, old_expressions))

  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X)) / 2)), max_depth)

  Z <- matrix(Z, nrow = length(X))
  XYZ <- cbind(scale(X), scale(Y), scale(Z))

  phi_x <- .condopt_marginal_likelihood(XYZ, target_idx = 1, z_idx = 3,
                                        z_min = 0, z_max = 1,
                                        c = 2 * c, rho = rho, depth = 1, max_depth, qdist)
  phi_y <- .condopt_marginal_likelihood(XYZ, target_idx = 2, z_idx = 3,
                                        z_min = 0, z_max = 1,
                                        c = 2 * c, rho = rho, depth = 1, max_depth, qdist)
  phi_xy_indep <- phi_x + phi_y

  phi_xy_dep <- .condopt_marginal_likelihood(XYZ, target_idx = c(1, 2), z_idx = 3,
                                             z_min = 0, z_max = 1,
                                             c = 1 * c, rho = rho, depth = 1, max_depth, qdist)

  bf <- phi_xy_indep - phi_xy_dep

  options(expressions = old_expressions)

  if (log_BF)
    return(list(bf = bf, p_H0 = NULL, p_H1 = NULL))

  bf <- exp(bf)
  return(list(bf = bf, p_H0 = 1 - 1 / (1 + bf), p_H1 = 1 / (1 + bf)))
}

.condopt_marginal_likelihood <- function(data, target_idx, z_idx, z_min, z_max,
                                         c, rho, depth, max_depth, qdist) {
  localData <- matrix(data[which(qdist(z_min) < data[, z_idx] & data[, z_idx] < qdist(z_max)),
                           target_idx], ncol = length(target_idx))
  logl <- .pt_marginal_likelihood(localData, low = rep(0, ncol(localData)),
                                         up = rep(1, ncol(localData)),
                                         c = c, depth = 1, max_depth, qdist)

  if (depth == max_depth || nrow(localData) <= 1) {
    return(logl)
  }

  phi_1 <- .condopt_marginal_likelihood(data, target_idx, z_idx, z_min, (z_min + z_max) / 2,
                                        c, rho, depth + 1, max_depth, qdist)
  phi_2 <- .condopt_marginal_likelihood(data, target_idx, z_idx, (z_min + z_max) / 2, z_max,
                                        c, rho, depth + 1, max_depth, qdist)

  return(matrixStats::logSumExp(c(log(rho) + logl, log(1 - rho) + phi_1 + phi_2)))
}

pt_independence_test <- function(X, Y, c = 1, max_depth = -1, qdist = qnorm, verbose = TRUE, log_BF = FALSE) {
  if (.is_discrete(X) || .is_discrete(Y)) {
    if (verbose)
      cat('Performing two-sample test\n')
    return(pt_d_sample_test(X, Y, c = c, max_depth = max_depth, qdist = qdist, log_BF))
  }

  if (verbose)
    cat('Performing continuous independence test\n')
  return(pt_continuous_independence_test(X, Y, c = c, max_depth = max_depth, qdist = qdist, log_BF))
}

# Implements a continuous independence test.
# See Filippi, S. and Holmes, C. C. (2017). A Bayesian 
# nonparametric approach to testing for dependence 
# between random variables. Bayesian Analysis, 12(4):919–938.
pt_continuous_independence_test <- function(X, Y, c = 1, max_depth = -1, qdist = qnorm, log_BF = FALSE) {
  old_expressions <- options()$expressions
  options(expressions = max(max_depth, old_expressions))

  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X)) / 2)), max_depth)

  X <- scale(X)
  Y <- scale(Y)
  XY <- cbind(X, Y)

  p_x <- .pt_marginal_likelihood(X, low = 0, up = 1, c = 2 * c,
                                        depth = 1, max_depth, qdist)
  p_y <- .pt_marginal_likelihood(Y, low = 0, up = 1, c = 2 * c,
                                        depth = 1, max_depth, qdist)
  p_xy <- .pt_marginal_likelihood(XY, low = c(0, 0), up = c(1, 1), c = 1 * c,
                                         depth = 1, max_depth, qdist)

  bf <- p_x + p_y - p_xy

  options(expressions = old_expressions)

  if (log_BF)
    return(list(bf = bf, p_H0 = NULL, p_H1 = NULL))

  bf <- exp(bf)
  return(list(bf = bf, p_H0 = 1 - 1 / (1 + bf), p_H1 = 1 / (1 + bf)))
}

# Implements a k-sample test, extension from:
# Holmes, C. C., Caron, F., Griffin, J. E., and Stephens, D. A. (2015).
# Two-sample Bayesian nonparametric hypothesis testing.
# Bayesian Analysis, 10(2):297–320.
pt_d_sample_test <- function(X, Y, c = 1, max_depth = -1, qdist = qnorm, log_BF = FALSE) {
  old_expressions <- options()$expressions
  options(expressions = max(max_depth, old_expressions))

  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X)) / 2)), max_depth)

  binary <- if (.is_discrete(X)) X else Y
  continuous <- if (.is_discrete(X)) Y else X

  data <- cbind(scale(continuous), binary)
  X <- data[, 1]

  p_H0 <- .pt_marginal_likelihood(X, low = 0, up = 1, c = c, depth = 1, max_depth, qdist)

  discrete_values <- if (length(unique(binary)) == 2) binary[1] else unique(binary)

  p_H1 <- max(sapply(discrete_values, function(i) {
    X1 <- data[data[, 2] == i, 1]
    X2 <- data[data[, 2] != i, 1]

    p_x1 <- .pt_marginal_likelihood(X1, low = 0, up = 1, c = c, depth = 1, max_depth, qdist)
    p_x2 <- .pt_marginal_likelihood(X2, low = 0, up = 1, c = c, depth = 1, max_depth, qdist)

    p_x1 + p_x2
  }))

  n_hypotheses <- length(discrete_values)
  bf <- p_H0 - p_H1 + log(n_hypotheses) # Bayes Factor with a Bonferroni-type correction

  options(expressions = old_expressions)

  if (log_BF)
    return(list(bf = bf, p_H0 = NULL, p_H1 = NULL))

  bf <- exp(bf)
  return(list(bf = bf, p_H0 = 1 - 1 / (1 + bf), p_H1 = 1 / (1 + bf)))
}

.pt_marginal_likelihood <- function(data, low, up, c, depth, max_depth, qdist) {
  if (depth == max_depth) {
    return(0)
  }

  if (length(low) == 1) {
    n_j <- c(length(which((qdist(low) < data) & (data <= qdist((low + up) / 2)))),
             length(which((qdist((low + up) / 2) < data) & (data <= qdist(up)))))
  } else {
    n_j <- c(length(which((qdist(low[1]) < data[, 1]) & (data[, 1] <= qdist((low[1] + up[1]) / 2)) &
                            (qdist(low[2]) < data[, 2]) & (data[, 2] <= qdist((low[2] + up[2]) / 2)))),
             length(which((qdist((low[1] + up[1]) / 2) < data[, 1]) & (data[, 1] <= qdist(up[1])) &
                            (qdist(low[2]) < data[, 2]) & (data[, 2] <= qdist((low[2] + up[2]) / 2)))),
             length(which((qdist(low[1]) < data[, 1]) & (data[, 1] <= qdist((low[1] + up[1]) / 2)) &
                            (qdist((low[2] + up[2]) / 2) < data[, 2]) & (data[, 2] <= qdist(up[2])))),
             length(which((qdist((low[1] + up[1]) / 2) < data[, 1]) & (data[, 1] <= qdist(up[1])) &
                            (qdist((low[2] + up[2]) / 2) < data[, 2]) & (data[, 2] <= qdist(up[2])))))
  }

  if (sum(n_j) == 0) {
    return(0)
  }

  a_j <- c * depth ^ 2

  if (length(n_j) == 2) {
    logl <- lbeta(n_j[1] + a_j, n_j[2] + a_j) - lbeta(a_j, a_j)
  } else {
    logl <- .lmbeta(n_j[1] + a_j, n_j[2] + a_j, n_j[3] + a_j, n_j[4] + a_j) -
      .lmbeta(a_j, a_j, a_j, a_j)
  }

  if (length(low) == 1) {
    likelihoods <- c(.pt_marginal_likelihood(data, low, (low + up) / 2, c, depth + 1, max_depth, qdist),
                     .pt_marginal_likelihood(data, (low + up) / 2, up, c, depth + 1, max_depth, qdist))
  } else {
    likelihoods <- c(.pt_marginal_likelihood(data, low, (low + up) / 2, c, depth + 1, max_depth, qdist),
                     .pt_marginal_likelihood(data, (low + up) / 2, up, c, depth + 1, max_depth, qdist),
                     .pt_marginal_likelihood(data,
                                                    c(low[1], (low[2] + up[2]) / 2),
                                                    c((low[1] + up[1]) / 2, up[2]),
                                                    c, depth + 1, max_depth, qdist),
                     .pt_marginal_likelihood(data,
                                                    c((low[1] + up[1]) / 2, low[2]),
                                                    c(up[1], (low[2] + up[2]) / 2),
                                                    c, depth + 1, max_depth, qdist))
  }

  return(logl + sum(likelihoods))
}

.lmbeta <- function(...) {
  sum(lgamma(c(...))) - lgamma(sum(c(...)))
}

.is_discrete <- function(X) {
  all(X %in% 0:10)
}
