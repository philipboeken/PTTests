# Maps
##############################################

### Interventions

ran <- c(2, 3, 4, 5, 6)

.mean_shift <- function(base, C) {
  theta <- sample(ran, 1)
  (1-C) * base + C * (base + theta)
}

.variance_shift <- function(base, C) {
  theta <- sample(ran, 1)
  (1-C) * base + C * (1+theta) * base
}

.fixed_point <- function(base, C) {
  theta <- sample(ran, 1)
  (1-C) * base + C * theta
}

.mixture <- function(base, C) {
  theta <- sample(ran, 1)
  idx <- sample(c(-1,theta), length(C), replace = TRUE)
  (1-C) * base + C * (base + idx)
}

.do_intervention <- function (int_options, base, C) {
  mapping <- sample(int_options, 1)[[1]]
  return(mapping(base, C))
}


### Link functions

.linear <- function(X) {
  X
}


.parabolic <- function(X) {
  X^2
}

.sinusoidal <- function(X) {
  sin(6*2*pi*X / (max(X)-min(X)))
}

.nonlin <- function(opts, X) {
  mapping <- sample(opts, 1)[[1]]
  return(mapping(X))
}

