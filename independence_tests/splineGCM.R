splineGCM <- function(x, y, S, data, verbose=FALSE) {
  suppressMessages(library(mgcv))
  
  if( verbose ) {
    cat( 'Testing ', x, '_||_', y, '|', S, '\n' )
  }
  stopifnot(length(x) == 1 && length(y) == 1)
  stopifnot(length(intersect(x,y)) == 0)
  stopifnot(length(intersect(union(x,y),S)) == 0)
  
  if( is.null(S) || length(S)==0 ) {
    if (all(data[,x] %in% 0:1)) {
      pval <- summary(gam(data[,x] ~ s(data[,y]), method = "REML", familiy='binomial'))$s.pv
    } else if (all(data[,y] %in% 0:1)) {
      pval <- summary(gam(data[,y] ~ s(data[,x]), method = "REML", familiy='binomial'))$s.pv
    } else {
      pval <- min(
        summary(gam(data[,y] ~ s(data[,x]), method = "REML"))$s.pv,
        summary(gam(data[,x] ~ s(data[,y]), method = "REML"))$s.pv
      )
    }
  } else {
    f <- sapply(S, function(col) s(data[,col]))
    f <- s(data[,S])
    if (all(data[,x] %in% 0:1)) {
      x_hat <- predict(gam(data[,x] ~ s(data[,S]), method = "REML", familiy='binomial'))
    } else {
      x_hat <- predict(gam(data[,x] ~ s(data[,S]), method = "REML"))
    }
    if (all(data[,y] %in% 0:1)) {
      y_hat <- predict(gam(data[,y] ~ s(data[,S]), method = "REML", familiy='binomial'))
    } else {
      y_hat <- predict(gam(data[,y] ~ s(data[,S]), method = "REML"))
    }
    pval <- gcm(x_hat, x, y_hat, y)$p.value
  }
  pval
}

gcm <- function(x_hat, x, y_hat, y) {
  R <- (x - x_hat) * (y - y_hat)
  T <- sqrt(length(x)) * mean(R) / sqrt(mean(R^2) - mean(R)^2)
  pval <- 2*pnorm(abs(T), lower.tail=FALSE)
  list(statistic = T, p.value = pval)
}
