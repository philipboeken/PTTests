source('independence_tests/bayesUCITest.R')
suppressWarnings(library(matrixStats))

bayes.CItest <- function (X, Y, Z = NULL, rho = 0.5, c = 1, max_depth = -1, qdist = qnorm, verbose = TRUE) {

  if (is.null(Z)) {
    if (verbose) {
      cat('Redirecting to UCI test\n')
    }
    return(bayes.UCItest(X, Y, c, max_depth, qdist, verbose))
  }

  max_depth <- ifelse(max_depth < 0, floor(log2(length(X))/2), max_depth)

  old_expressions <- options()$expressions
  options(expressions=2*max(max_depth, old_expressions))

  if (all(X %in% 0:1) && all(Y %in% 0:1) || length(X) <= 2) {
    return(list(bf=1, p_H0=1/2, p_H1=1/2))
  }
  
  if (all(X %in% 0:1) || all(Y %in% 0:1)) {

    if (verbose) {
      cat('Performing two-sample test\n')
    }
    
    if (all(X %in% 0:1)) {
      bin <- X
      cont <- Y
    } else {
      cont <- X
      bin <- Y
    }

    data <- cbind(scale(cont), bin, scale(Z))
    XZ <- data[, c(1, 3)]
    X1Z <- data[data[,2] == 0, c(1, 3)]
    X2Z <- data[data[,2] == 1, c(1, 3)]
    
    p_x <- .opt_pt(XZ, 1, 2, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    p_x1 <- .opt_pt(X1Z, 1, 2, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    p_x2 <- .opt_pt(X2Z, 1, 2, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)

    bf <- exp(p_x - p_x1 - p_x2)

  } else {

    XYZ <- cbind(scale(X), scale(Y), scale(Z))

    phi_x <- .opt_pt(XYZ, 1, 3, z_min=0, z_max=1, c, rho, depth=1, max_depth, qdist)
    phi_y <- .opt_pt(XYZ, 2, 3, z_min=0, z_max=1, c=c, rho=rho, depth=1, max_depth, qdist)
    phi_xy <- .opt_pt(XYZ, c(1, 2), 3, z_min=0, z_max=1, c=c, rho=rho, depth=1, max_depth, qdist)

    bf <- exp(phi_x + phi_y - phi_xy)

  }

  options(expressions = old_expressions)

  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}


.opt_pt <- function(data, target_idx, z_idx, z_min, z_max, c, rho, depth, max_depth, qdist) {
  local <- matrix(data[which(qdist(z_min) < data[, z_idx] & data[, z_idx] < qdist(z_max)),
                       target_idx], ncol=length(target_idx))
  logl <- .pt_logl(local, low=rep(0, ncol(local)), up=rep(1, ncol(local)),
                  c=c, depth=1, max_depth, qdist)

  if (depth == max_depth || nrow(local) <= 1) {
    return(logl)
  }

  phi_1 <- .opt_pt(data, target_idx, z_idx, z_min, (z_min + z_max)/2,
                  c, rho, depth + 1, max_depth, qdist)
  phi_2 <- .opt_pt(data, target_idx, z_idx, (z_min + z_max)/2, z_max,
                  c, rho, depth + 1, max_depth, qdist)
  
  return(logSumExp(c(log(rho) + logl, log(1-rho) + phi_1 + phi_2)))
}
