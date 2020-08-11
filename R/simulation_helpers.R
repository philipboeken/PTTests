# Maps
##############################################

### Interventions

no_intervention <- function(base, C) {
  base
}

mean_shift <- function(base, C) {
  base + C
}

variance_shift <- function(base, C) {
  base * (1 + C)
}

fixed_point <- function(base, C) {
  C
}

mixture <- function(base, C) {
  base + (C > max(C)/2)
}

.do_intervention <- function (int_options, base, C) {
  mapping <- sample(int_options, 1)[[1]]
  return(mapping(base, C))
}


### Link functions

linear <- function(X) {
  X
}


parabolic <- function(X) {
  X^2
}

sinusoidal <- function(X) {
  sin(2*pi*X / sd(X))
}

.nonlin <- function(opts, X) {
  mapping <- sample(opts, 1)[[1]]
  return(mapping(X))
}

### Generate data ##########################

.graph_1 <- function (n, dim_C, err_sd, p_link, interv_options, nonlin_options) {
  # C -> X -> Y
  
  C <- sample(1:dim_C, n, replace = TRUE)
  
  intervene <- rbinom(1, 1, p_link)
  X <- intervene * .do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
  
  link_nonlin <- rbinom(1, 1, p_link)
  Y <- link_nonlin * .nonlin(nonlin_options, X)
  Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
  
  list(C = C, X = X, Y = Y, label_CX = 1 - intervene, label_XY = 1 - link_nonlin, 
       label_CY_X = 1, label_lcd = intervene * link_nonlin)
}

.graph_2 <- function (n, dim_C, err_sd, p_link, interv_options, nonlin_options) {
  # C -> X <- Y
  
  C <- sample(1:dim_C, n, replace = TRUE)
  
  Y <- rnorm(n)
  
  link_nonlin <- rbinom(1, 1, p_link)
  X <- link_nonlin * .nonlin(nonlin_options, Y)
  X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
  
  intervene <- rbinom(1, 1, p_link)
  X <- intervene * .do_intervention(interv_options, X, C) + (1-intervene) * X
  
  list(C = C, X = X, Y = Y, label_CX = 1 - intervene, label_XY = 1 - link_nonlin, 
       label_CY_X = 1 - link_nonlin * intervene, label_lcd = 0)
}

.graph_3 <- function (n, dim_C, err_sd, p_link, interv_options, nonlin_options) {
  # C -> X <- L -> Y
  
  C <- sample(1:dim_C, n, replace = TRUE)
  
  L <- rnorm(n)
  
  link_nonlin1 <- rbinom(1, 1, p_link)
  Y <- link_nonlin1 * .nonlin(nonlin_options, L)
  Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
  
  link_nonlin2 <- rbinom(1, 1, p_link)
  X <- link_nonlin2 * .nonlin(nonlin_options, L)
  X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
  
  intervene <- rbinom(1, 1, p_link)
  X <- intervene * .do_intervention(interv_options, X, C) + (1-intervene) * X
  
  link_nonlin <- link_nonlin1 & link_nonlin2
  
  list(C = C, X = X, Y = Y, label_CX = 1 - intervene, label_XY = 1 - link_nonlin, 
       label_CY_X = 1 - link_nonlin * intervene, label_lcd = 0)
}

get_data <- function(graph_probs, n, dim_C, err_sd, p_link, interv_options, nonlin_options) {
  graph <- sample(c(.graph_1, .graph_2, .graph_3), 1, prob=graph_probs)[[1]]
  
  graph(n, dim_C, err_sd, p_link, interv_options, nonlin_options)
}

