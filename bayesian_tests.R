source('indepTests/bayesCITest.R')

data_1 <- function(n) {
  return(list(X=runif(n), Y=runif(n), Z=runif(n)))
}

data_2 <- function(n) {
  U <- runif(n)
  V <- runif(n)
  Z <- runif(n)
  return(list(X=(U+Z)/2, Y=(V+Z)/2, Z=Z))
}

data_3 <- function(n) {
  X <- runif(n)
  Y <- runif(n)
  W <- rnorm(n, 0, 1/2)
  Z <- atan(Y / X + W) / pi + 1/2
  return(list(X=X, Y=Y, Z=Z))
}

data_4 <- function(n) {
  suppressWarnings(library(mvtnorm))
  mv <- rmvnorm(n, mean=c(0.5, 0.5), sigma=cbind(c(2, 1), c(1, 2))/20)
  mv[((mv[,1] < 0) | (1 < mv[,1])) & ((mv[,2] < 0) | (1 < mv[,2])), ] <- 0
  X <- mv[,1]
  Y <- mv[,2]
  B <- rbinom(n, 1, 0.9)
  Z <- rnorm(n, 1/2, 1/2)
  Z[(Z < 0) | (1 < Z)] <- 0
  Z <- Z * (B + (1-B) * X * Y)
  return(list(X=X, Y=Y, Z=Z))
}

n <- 500

d <- data_3(n)
X <- d$X
Y <- d$Y
Z <- d$Z
# result <- bayes.UCItest(X, Y, max_dept=-1)
# cat("\nBF(H0, H1):\t", result$bf, "\nP(H1 | XYZ):\t", result$p_H1, "\n")
# result <- bayes.CItest(X, Y, Z, rho=0.5)
# cat("\nBF(H0, H1):\t", result$bf, "\nP(H1 | XYZ):\t", result$p_H1, "\n")

# res <- c()
# for (i in 1:100) {
#   # C <- c(rep(0,n), rep(1, n))
  C <- rbinom(n, 1, 0.5)
  X <- (1 - C)*rnorm(n)
  result <- bayes.UCItest(X, C, max_depth = -1, verbose=FALSE)
  cat(result$p_H1)
  # res <- c(res, result$p_H1)
# }
# hist(res)

