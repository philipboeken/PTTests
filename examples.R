source('init.R', chdir=TRUE)

esd <- 0.5
alpha <- 0.05
n <- 400

# nonlin <- function(x) sin(3*x)
nonlin <- function(x) x^2

results <- data.frame(
  test=character(0),
  true=integer(0),
  linear=integer(0),
  pcor=numeric(0),
  splineGCM=numeric(0),
  RCoT=numeric(0),
  bayes=numeric(0))

add_results <- function(dgp, true, linear, pcor, splineGCM, RCoT, bayes) {
  rbind(results, data.frame(dgp=dgp, true=true, linear=linear,
    pcor=round(pcor, digits=3), splineGCM=round(splineGCM, digits=3), RCoT=round(RCoT, digits=3),
    bayes=round(bayes, digits=3)))
}

dgp <- "X <- S2 <- S1 -> Y (X_||_Y|{S1,S2})"
# LINEAR
S1 <- rnorm(n)
S2 <- S1 + rnorm(n, sd=esd*sd(S1))
X <- S2 + rnorm(n, sd=esd*sd(S2))
Y <- S1 + rnorm(n, sd=esd*sd(S1))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)
# NONLINEAR
S1 <- rnorm(n)
S2 <- S1 + rnorm(n, sd=esd*sd(S1))
X <- nonlin(S2) + rnorm(n, sd=esd*sd(nonlin(S2)))
Y <- nonlin(S1) + rnorm(n, sd=esd*sd(nonlin(S1)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)

dgp <- "X -> S2 <- S1 -> Y (X_||_Y|{S1,S2})"
# LINEAR
X <- rnorm(n)
S1 <- rnorm(n)
Y <- S1 + rnorm(n, sd=esd*sd(S1))
S2 <- S1 + X + rnorm(n,sd=esd*sd(S1+X))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)
# NONLINEAR, X~S2 should not work
X <- rnorm(n)
S1 <- rnorm(n)
Y <- nonlin(S1) + rnorm(n, sd=esd*sd(nonlin(S1)))
S2 <- S1 + nonlin(X) + rnorm(n,sd=esd*sd(S1 + nonlin(X)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)

dgp <- "X -> S2 <- S1 <- Y (X_||_Y|{S1,S2})"
# LINEAR
X <- rnorm(n)
Y <- rnorm(n)
S1 <- Y + rnorm(n, sd=esd*sd(Y))
S2 <- S1 + X + rnorm(n,sd=esd*sd(S1+X))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)
# NONLINEAR, X~S2 should not work
X <- rnorm(n)
Y <- rnorm(n)
S1 <- nonlin(Y) + rnorm(n, sd=esd*sd(nonlin(Y)))
S2 <- S1 + nonlin(X) + rnorm(n,sd=esd*sd(S1 + nonlin(X)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)

dgp <- "X -> S2 -> S1 -> Y (X_||_Y|{S1,S2})"
# LINEAR
X <- rnorm(n)
S2 <- X + rnorm(n,sd=esd*sd(X))
S1 <- S2 + rnorm(n, sd=esd*sd(S2))
Y <- S1 + rnorm(n, sd=esd*sd(S1))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)
# NONLINEAR, X~S2 should not work
X <- rnorm(n)
S2 <- nonlin(X) + rnorm(n,sd=esd*sd(nonlin(X)))
S1 <- nonlin(S2) + rnorm(n,sd=esd*sd(nonlin(S2)))
Y <- nonlin(S1) + rnorm(n,sd=esd*sd(nonlin(S1)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, cbind(S1, S2))$p,
                       splineGCM(1, 2, c(3, 4), cbind(X, Y, S1, S2)),
                       RCoT(X, Y, cbind(S1, S2))$p,
                       NaN)

dgp <- "X <- Z -> Y (X_||_Y|Z)"
# LINEAR
Z <- rnorm(n)
X <- Z + rnorm(n, sd=esd*sd(Z))
Y <- Z + rnorm(n, sd=esd*sd(Z))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)
# NONLINEAR
Z <- rnorm(n)
X <- nonlin(Z) + rnorm(n, sd=esd*sd(nonlin(Z)))
Y <- nonlin(Z) + rnorm(n, sd=esd*sd(nonlin(Z)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)

dgp <- "X -> Z <- Y (not X_||_Y|Z)"
# LINEAR
X <- rnorm(n)
Y <- rnorm(n)
Z <- X + Y + rnorm(n, sd=esd*sd(X+Y))
results <- add_results(dgp, 0, 'linear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)
# NONLINEAR
X <- rnorm(n)
Y <- rnorm(n)
Z <- nonlin(X) + nonlin(Y) + rnorm(n, sd=esd*sd(nonlin(X) + nonlin(Y)))
pcor.test(X, Y, Z)$p
results <- add_results(dgp, 0, 'nonlinear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)

dgp <- "X -> Z -> Y (X_||_Y|Z)"
# LINEAR
X <- rnorm(n)
Z <- X + rnorm(n, sd=esd*sd(X))
Y <- Z + rnorm(n, sd=esd*sd(Z))
results <- add_results(dgp, 1, 'linear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)
# NONLINEAR, X~Z should not work
X <- rnorm(500)
Z <- nonlin(X) + rnorm(500, sd=esd*sd(nonlin(X)))
Y <- nonlin(Z) + rnorm(500, sd=esd*sd(nonlin(Z)))
results <- add_results(dgp, 1, 'nonlinear',
                       pcor.test(X, Y, Z)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z)),
                       RCoT(X, Y, Z)$p,
                       bayes.CItest(X, Y, Z)$p_H0)

dgp <- "Z2 -> X -> Z1 <- Y <- Z2 (not X_||_Y|Z1)"
# LINEAR
Z2 <- rnorm(n)
X <- 0.2*Z2 + rnorm(n, sd=esd*sd(0.2*Z2))
Y <- 0.2*Z2 + rnorm(n, sd=esd*sd(0.2*Z2))
Z1 <- X + Y + rnorm(n, sd=esd*sd(X + Y))
results <- add_results(dgp, 0, 'linear',
                       pcor.test(X, Y, Z1)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z1)),
                       RCoT(X, Y, Z1)$p,
                       bayes.CItest(X, Y, Z1)$p_H0)
# NONLINEAR
Z2 <- rnorm(n)
X <- nonlin(Z2) + rnorm(n, sd=esd*sd(nonlin(Z2)))
Y <- nonlin(Z2) + rnorm(n, sd=esd*sd(nonlin(Z2)))
Z1 <- X + Y + rnorm(n, sd=esd*sd(X + Y))
results <- add_results(dgp, 0, 'nonlinear',
                       pcor.test(X, Y, Z1)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z1)),
                       RCoT(X, Y, Z1)$p,
                       bayes.CItest(X, Y, Z1)$p_H0)

dgp <- "Z2 -> X <- Z1 -> Y <- Z2 (not X_||_Y|Z1)"
# LINEAR
Z1 <- rnorm(n)
Z2 <- rnorm(n)
X <- Z1 + 0.2*Z2 + rnorm(n, sd=esd*sd(Z1 + 0.2*Z2))
Y <- Z1 + 0.2*Z2 + rnorm(n, sd=esd*sd(Z1 + 0.2*Z2))
results <- add_results(dgp, 0, 'linear',
                       pcor.test(X, Y, Z1)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z1)),
                       RCoT(X, Y, Z1)$p,
                       bayes.CItest(X, Y, Z1)$p_H0)
# NONLINEAR
Z1 <- rnorm(n)
Z2 <- rnorm(n)
X <- nonlin(Z1) + 0.2*Z2 + rnorm(n, sd=esd*sd(nonlin(Z1) + 0.2*Z2))
Y <- nonlin(Z1) + 0.2*Z2 + rnorm(n, sd=esd*sd(nonlin(Z1) + 0.2*Z2))
results <- add_results(dgp, 0, 'nonlinear',
                       pcor.test(X, Y, Z1)$p,
                       splineGCM(1, 2, c(3), cbind(X, Y, Z1)),
                       RCoT(X, Y, Z1)$p,
                       bayes.CItest(X, Y, Z1)$p_H0)

# C -> X -> Y, LCD
C <- rbinom(n, 1, 0.5)
C <- c(rep(0,n), rep(1, n))
X <- 2*C + rexp(2*n)
Y <- nonlin(X) + rnorm(2*n, sd=esd*sd(nonlin(X)))
# C_||_X
cor.test(C, X)$p.value
splineGCM(1, 2, c(), cbind(C, X))
RCoT(C, X)$p
bayes.UCItest(X, C, max_depth = -1)$p_H0
# X_||_Y
cor.test(X, Y)$p.value
splineGCM(1, 2, c(), cbind(X, Y))
RCoT(X, Y)$p
bayes.UCItest(X, Y, max_depth=-1)$p_H0
# C_||_Y|X
splineGCM(1, 2, c(3), cbind(C, Y, X))
pcor.test(C, Y, X)$p
RCoT(C, Y, X)$p
bayes.CItest(Y, C, X, max_depth=-1)$p_H0
