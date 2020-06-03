polyatree.CITest <- function (X, Y, Z=NULL, rho=0.5, c=1, max_depth=-1, qdist=qnorm, verbose=TRUE) {
  
  if (is.null(Z)) {
    if (verbose)
      cat('Redirecting to independence test\n')
    return(polyatree.IndepTest(X, Y, c, max_depth, qdist, verbose))
  }
  
  if (all(X %in% 0:1) || all(Y %in% 0:1)) {
    if (verbose)
      cat('Performing conditional two-sample test\n')
    return(polyatree.TwoSampleCITest(X, Y, Z, rho=rho, c=c, max_depth=max_depth, qdist=qdist))
    
  }
  
  if (verbose)
    cat('Performing continuous conditional independence test\n')
  return(polyatree.ContinuousCITest(X, Y, Z, rho=rho, c=c, max_depth=max_depth, qdist=qdist))
}

polyatree.TwoSampleCITest <- function (X, Y, Z, rho=0.5, c=1, max_depth=-1, qdist=qnorm) {
  old_expressions <- options()$expressions
  options(expressions=max(max_depth, old_expressions))
  
  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X))/2)), max_depth)
  
  if (all(X %in% 0:1) && all(Y %in% 0:1) || length(X) <= 2) {
    return(list(bf=1, p_H0=1/2, p_H1=1/2))
  }
  
  if (all(X %in% 0:1)) {
    binary <- X
    continuous <- Y
  } else {
    binary <- Y
    continuous <- X
  }
  data <- cbind(scale(continuous), binary, scale(Z))
  XZ <- data[, c(1, 3)]
  X1Z <- data[data[,2] == 0, c(1, 3)]
  X2Z <- data[data[,2] == 1, c(1, 3)]
  
  p_x <- .condOPT.MarginalLikelihood(XZ, 1, 2, z_min=0, z_max=1, c=c, 
                                     rho=rho, depth=1, max_depth, qdist)
  p_x1 <- .condOPT.MarginalLikelihood(X1Z, 1, 2, z_min=0, z_max=1, c=c, 
                                      rho=rho, depth=1, max_depth, qdist)
  p_x2 <- .condOPT.MarginalLikelihood(X2Z, 1, 2, z_min=0, z_max=1, c=c, 
                                      rho=rho, depth=1, max_depth, qdist)
  
  bf <- exp(p_x - p_x1 - p_x2)
  
  options(expressions=old_expressions)
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}

polyatree.ContinuousCITest <- function (X, Y, Z=NULL, rho=0.5, c=1, 
                                        max_depth=-1, qdist=qnorm, verbose=TRUE) {
  old_expressions <- options()$expressions
  options(expressions=max(max_depth, old_expressions))
  
  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X))/2)), max_depth)
  
  XYZ <- cbind(scale(X), scale(Y), scale(Z))
  
  phi_x <- .condOPT.MarginalLikelihood(XYZ, 1, 3, z_min=0, z_max=1, c=2*c, 
                                       rho=rho, depth=1, max_depth, qdist)
  phi_y <- .condOPT.MarginalLikelihood(XYZ, 2, 3, z_min=0, z_max=1, c=2*c, 
                                       rho=rho, depth=1, max_depth, qdist)
  phi_xy <- .condOPT.MarginalLikelihood(XYZ, c(1, 2), 3, z_min=0, z_max=1, 
                                        c=1*c, rho=rho, depth=1, max_depth, qdist)
  
  bf <- exp(phi_x + phi_y - phi_xy)
  
  options(expressions=old_expressions)
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}

.condOPT.MarginalLikelihood <- function(data, target_idx, z_idx, z_min, z_max,
                                        c, rho, depth, max_depth, qdist) {
  localData <- matrix(data[which(qdist(z_min) < data[, z_idx] & data[, z_idx] < qdist(z_max)),
                           target_idx], ncol=length(target_idx))
  logl <- .pt.MarginalLikelihood(localData, low=rep(0, ncol(localData)), up=rep(1, ncol(localData)),
                                 c=c, depth=1, max_depth, qdist)
  
  if (depth == max_depth || nrow(localData) <= 1) {
    return(logl)
  }
  
  phi_1 <- .condOPT.MarginalLikelihood(data, target_idx, z_idx, z_min, (z_min + z_max)/2,
                                       c, rho, depth + 1, max_depth, qdist)
  phi_2 <- .condOPT.MarginalLikelihood(data, target_idx, z_idx, (z_min + z_max)/2, z_max,
                                       c, rho, depth + 1, max_depth, qdist)
  
  return(matrixStats::logSumExp(c(log(rho) + logl, log(1-rho) + phi_1 + phi_2)))
}

polyatree.IndepTest <- function (X, Y, c=1, max_depth=-1, qdist=qnorm, verbose=TRUE) {
  if (all(X %in% 0:1) || all(Y %in% 0:1)) {
    if (verbose)
      cat('Performing two-sample test\n')
    return(polyatree.TwoSampleTest(X, Y, c =c, max_depth=max_depth, qdist=qdist))
  }
  
  if (verbose)
    cat('Performing continuous independence test\n')
  return(polyatree.ContinuousIndepTest(X, Y, c =c, max_depth=max_depth, qdist=qdist))
}

polyatree.ContinuousIndepTest <- function (X, Y, c=1, max_depth=-1, qdist=qnorm) {
  old_expressions <- options()$expressions
  options(expressions=max(max_depth, old_expressions))
  
  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X))/2)), max_depth)
  
  X <- scale(X)
  Y <- scale(Y)
  XY <- cbind(X, Y)
  
  p_x <- .pt.MarginalLikelihood(X, low=0, up=1, c=2*c, depth=1, max_depth, qdist)
  p_y <- .pt.MarginalLikelihood(Y, low=0, up=1, c=2*c, depth=1, max_depth, qdist)
  p_xy <- .pt.MarginalLikelihood(XY, low=c(0, 0), up=c(1, 1), c=1*c, depth=1, max_depth, qdist)
  
  bf <- exp(p_x + p_y - p_xy)
  
  options(expressions=old_expressions)
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}

polyatree.TwoSampleTest <- function(X, Y, c=1, max_depth=-1, qdist=qnorm) {
  old_expressions <- options()$expressions
  options(expressions=max(max_depth, old_expressions))
  
  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X))/2)), max_depth)
  
  if (all(X %in% 0:1)) {
    binary <- X
    continuous <- Y
  } else {
    binary <- Y
    continuous <- X
  }
  data <- cbind(scale(continuous), binary)
  X <- data[, 1]
  X1 <- data[data[,2] == 0, 1]
  X2 <- data[data[,2] == 1, 1]
  
  p_xy <- .pt.MarginalLikelihood(X, low=0, up=1, c=c, depth=1, max_depth, qdist)
  p_x <- .pt.MarginalLikelihood(X1, low=0, up=1, c=c, depth=1, max_depth, qdist)
  p_y <- .pt.MarginalLikelihood(X2, low=0, up=1, c=c, depth=1, max_depth, qdist)
  
  bf <- exp(p_xy - p_x - p_y)
  
  options(expressions=old_expressions)
  
  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}


.pt.MarginalLikelihood <- function(data, low, up, c, depth, max_depth, qdist) {
  if (depth == max_depth) {
    return(0)
  }
  
  if (length(low) == 1) {
    n_j <- c(length(which((qdist(low) < data) & (data < qdist((low + up)/2)))),
             length(which((qdist((low + up)/2) < data) & (data < qdist(up)))))
  } else {
    n_j <- c(length(which((qdist(low[1]) < data[,1]) & (data[,1] < qdist((low[1] + up[1])/2)) &
                            (qdist(low[2]) < data[,2]) & (data[,2] < qdist((low[2] + up[2])/2)))),
             length(which((qdist((low[1] + up[1])/2) < data[,1]) & (data[,1] < qdist(up[1])) &
                            (qdist(low[2]) < data[,2]) & (data[,2] < qdist((low[2] + up[2])/2)))),
             length(which((qdist(low[1]) < data[,1]) & (data[,1] < qdist((low[1] + up[1])/2)) &
                            (qdist((low[2] + up[2])/2) < data[,2]) & (data[,2] < qdist(up[2])))),
             length(which((qdist((low[1] + up[1])/2) < data[,1]) & (data[,1] < qdist(up[1])) &
                            (qdist((low[2] + up[2])/2) < data[,2]) & (data[,2] < qdist(up[2])))))
  }
  
  if (sum(n_j) < 1) {
    return(0)
  }
  
  a_j <- c * depth^2
  if (length(n_j) == 2) {
    logl <- lbeta(n_j[1] + a_j, n_j[2] + a_j) - lbeta(a_j, a_j)
  } else {
    logl <- .lmbeta(n_j[1] + a_j, n_j[2] + a_j, n_j[3] + a_j, n_j[4] + a_j) -
      .lmbeta(a_j, a_j, a_j, a_j)
  }
  
  if (length(low) == 1) {
    likelihoods <- c(.pt.MarginalLikelihood(data, low, (low + up)/2, c, depth + 1, max_depth, qdist),
                     .pt.MarginalLikelihood(data, (low + up)/2, up, c, depth + 1, max_depth, qdist))
  } else {
    likelihoods <- c(.pt.MarginalLikelihood(data, low, (low + up)/2, c, depth + 1, max_depth, qdist),
                     .pt.MarginalLikelihood(data, (low + up)/2, up, c, depth + 1, max_depth, qdist),
                     .pt.MarginalLikelihood(data, 
                                            c(low[1], (low[2] + up[2])/2), 
                                            c((low[1] + up[1])/2, up[2]),
                                            c, depth + 1, max_depth, qdist),
                     .pt.MarginalLikelihood(data, 
                                            c((low[1] + up[1])/2, low[2]), 
                                            c(up[1], (low[2] + up[2])/2),
                                            c, depth + 1, max_depth, qdist))
  }
  
  return(logl + sum(likelihoods))
}

.lmbeta <- function(...) {
  return(sum(lgamma(c(...)))-lgamma(sum(c(...))))
}
