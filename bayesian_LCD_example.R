source('indepTests/bayesCITest.R')

# LCD EXAMPLE WHERE ALL TESTS WORK PROPERLY
nonlin <- function(x) x^2
nonlin <- function(x) sin(3*x)
n <- 200
esd <- 0.1

# C <- c(rep(0,n/2), rep(1, n/2))
C <- rbinom(n, 1, 0.5)
X <- (1+2*C) * rexp(n)

result <- bayes.UCItest(X, C, max_depth = -1)
cat("\nBF(H0, H1):\t", result$bf, "\nP(H0 | XY):\t", result$p_H0, "\n")

Y <- nonlin(X) + rnorm(n, sd=esd*sd(nonlin(X)))
# Y <- 3*C + rnorm(n, sd=esd)
result <- bayes.UCItest(X, Y, max_depth = -1)
cat("\nBF(H0, H1):\t", result$bf, "\nP(H0 | XY):\t", result$p_H0, "\n")

result <- bayes.CItest(Y, C, X, max_depth = -1)
cat("\nBF(H0, H1):\t", result$bf, "\nP(H0 | Z):\t", result$p_H0, "\n")

# Example of output:
# Performing two-sample test
# BF(H0, H1):	 1.171266e-69 
# P(H0 | XY):	 0 
# 
# BF(H0, H1):	 1.0208e-09 
# P(H0 | XY):	 1.0208e-09 
#
# Performing two-sample test
# BF(H0, H1):	 624086.4 
# P(H0 | Z):	 0.9999984 
