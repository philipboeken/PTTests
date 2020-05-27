bayes.UCItest <- function (X, Y, c = 1, max_depth = -1, qdist = qnorm, verbose=TRUE) {

  max_depth <- ifelse(max_depth < 0, max(1, floor(log2(length(X))/2)), max_depth)

  old_expressions <- options()$expressions
  options(expressions = max(max_depth, old_expressions))

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

    data <- cbind(scale(cont), bin)
    X <- data[, 1]
    X1 <- data[data[,2] == 0, 1]
    X2 <- data[data[,2] == 1, 1]

    p_xy <- .pt_logl(X, low=0, up=1, c=1*c, depth=1, max_depth, qdist)
    p_x <- .pt_logl(X1, low=0, up=1, c=1*c, depth=1, max_depth, qdist)
    p_y <- .pt_logl(X2, low=0, up=1, c=1*c, depth=1, max_depth, qdist)

    bf <- exp(p_xy - p_x - p_y)
    
  } else {
    
    X <- scale(X)
    Y <- scale(Y)
    XY <- cbind(X, Y)

    p_x <- .pt_logl(X, low=0, up=1, c=2*c, depth=1, max_depth, qdist)
    p_y <- .pt_logl(Y, low=0, up=1, c=2*c, depth=1, max_depth, qdist)
    p_xy <- .pt_logl(XY, low=c(0, 0), up=c(1, 1), c=1*c, depth=1, max_depth, qdist)

    bf <- exp(p_x + p_y - p_xy)
    
  }

  options(expressions = old_expressions)

  return(list(bf=bf, p_H0=1-1/(1+bf), p_H1=1/(1+bf)))
}


.pt_logl <- function(data, low, up, c, depth, max_depth, qdist) {
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

  logl <- .logl_j(n_j, a_j=c * depth^2)

  if (length(low) == 1) {
    sub_pt_logl <- c(.pt_logl(data, low, (low + up)/2, c, depth + 1, max_depth, qdist),
                     .pt_logl(data, (low + up)/2, up, c, depth + 1, max_depth, qdist))
  } else {
    sub_pt_logl <- c(.pt_logl(data, low, (low + up)/2, c, depth + 1, max_depth, qdist),
                     .pt_logl(data, (low + up)/2, up, c, depth + 1, max_depth, qdist),
                     .pt_logl(data, c(low[1], (low[2] + up[2])/2), c((low[1] + up[1])/2, up[2]),
                             c, depth + 1, max_depth, qdist),
                     .pt_logl(data, c((low[1] + up[1])/2, low[2]), c(up[1], (low[2] + up[2])/2),
                             c, depth + 1, max_depth, qdist))
  }

  return(logl + sum(sub_pt_logl))
}

.logl_j <- function(n_j, a_j) {
  if (length(n_j) == 2) {
    return(lbeta(n_j[1] + a_j, n_j[2] + a_j) - lbeta(a_j, a_j))
  } else {
    return(.lmbeta(n_j[1] + a_j, n_j[2] + a_j, n_j[3] + a_j, n_j[4] + a_j) -
             .lmbeta(a_j, a_j, a_j, a_j))
  }
}

.lmbeta <- function(...) {
  return(sum(lgamma(c(...)))-lgamma(sum(c(...))))
}
