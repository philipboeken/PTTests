# Clear workspace
rm(list = ls(all.names = TRUE))

# Imports
source('independence_tests/test_wrappers.R')
source('independence_tests/RhoBFP.R')
source('experiments/simulations/maps.R')
source('helpers.R')
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))
suppressWarnings(library(cowplot))
library(latex2exp)


# Input parameters
##############################################
n <- 300
m <- 500

err_sd <- 0.1

p_ci <- 0.6
p_link <- 0.8
p_two_sample <- 0.5

nonlin_options <- c(
  linear,
  parabolic,
  sinusoidal
)

interv_options <- c(
  mean_shift,
  variance_shift,
  fixed_point,
  mixture
)


# Setup test
##############################################
get_data <- function(n, p_two_sample, p_link, p_ci, err_sd, nonlin_options, interv_options) {
  C <- rbinom(n, 1, p_two_sample)
  
  cond_indep <- rbinom(1, 1, p_ci)
  if (cond_indep) { # C -> Z -> X
    intervene <- rbinom(1, 1, p_link)
    Z <- intervene * do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
    
    link_nonlin <- rbinom(1, 1, p_link)
    X <- link_nonlin * nonlin(nonlin_options, Z)
    X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
  } else {
    if (runif(1) <= 0.5) { # C -> Z <- X
      X <- rnorm(n)
      
      link_nonlin <- rbinom(1, 1, p_link)
      Z <- link_nonlin * nonlin(nonlin_options, X)
      Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      
      intervene <- rbinom(1, 1, p_link)
      Z <- intervene * do_intervention(interv_options, Z, C) + (1-intervene) * Z
    } else { # C -> Z <- L -> X
      L <- rnorm(n)
      
      link_nonlin1 <- rbinom(1, 1, p_link)
      X <- link_nonlin1 * nonlin(nonlin_options, L)
      X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
      
      link_nonlin2 <- rbinom(1, 1, p_link)
      Z <- link_nonlin2 * nonlin(nonlin_options, L)
      Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      
      intervene <- rbinom(1, 1, p_link)
      Z <- intervene * do_intervention(interv_options, Z, C) + (1-intervene) * Z
      
      link_nonlin <- link_nonlin1 & link_nonlin2
    }
  }
  
  cond_indep <- as.numeric(cond_indep | !link_nonlin | !intervene)
  lcd <- as.numeric(intervene & link_nonlin & cond_indep)
  
  return(list(C=C, X=X, Z=Z,
              label_ts=as.numeric(intervene),
              label_uci=as.numeric(link_nonlin),
              label_ci=cond_indep,
              label_lcd=lcd))
}

get_results <- function(dataset, test, bayesian=FALSE){
  result <- foreach(i=1:length(dataset), .combine=rbind) %dopar% {
    data <- dataset[[i]]
    ts <- test(data$C, data$Z)
    uci <- test(data$Z, data$X)
    ci <- test(data$C, data$X, data$Z)
    CZ <- test(data$C, data$Z)
    
    if (bayesian) {
      lcd_CY <- (ts <= 1/2) * (uci <= 1/2) * (ci > 1/2) * (1-CZ)
    } else {
      n <- length(data$C)
      lcd_CY <- (ts <= 1/(5*sqrt(n))) * (uci <= 1/(5*sqrt(n))) * (ci > 1/(5*sqrt(n))) * (1-CZ)
    }
    
    lcd_min <- min((1 - ts), (1 - uci), ci)
    
    return(data.frame(
      label_ts=data$label_ts,
      label_uci=data$label_uci,
      label_ci=data$label_ci,
      label_CY=data$label_ci,
      label_lcd_min=data$label_lcd,
      label_lcd_CY=data$label_lcd,
      ts=1-ts,
      uci=1-uci,
      ci=ci,
      CY=1-CZ,
      lcd_min=lcd_min,
      lcd_CY=lcd_CY
    ))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)
data <- lapply(1:m, function (i) get_data(n, p_two_sample, p_link, p_ci, 
                                          err_sd, nonlin_options, interv_options))
results <- list(
  # pcor=get_results(data, .pcor_wrapper),
  bayes=get_results(data, .bayes_wrapper, TRUE),
  # bcor_pb=get_results(data, .bayes_transform(.pcor_wrapper))
  # bcor=get_results(data, .bcor_wg_wrapper, TRUE),
  # bcor_approx=get_results(data, .bcor_approx_wrapper),
  # bcor_ly=get_results(data, .bcor_ly_wrapper),
  # gcm_bayes=get_results(data, .bayes_transform(.gcm_wrapper)),
  gcm=get_results(data, .gcm_wrapper)
  # ccit=get_results(data, .ccit_wrapper),
  # ccit_bayes=get_results(data, .bayes_transform(.ccit_wrapper)),
  # rcot=get_results(data, .rcot_wrapper)
  # rcot_bayes=get_results(data, .bayes_transform(.rcot_wrapper))
)

stopCluster(.cl)


# Process results
##############################################

.get_results_by_type <- function (results, type) {
  result <- data.frame(label=results$bayes[,{{paste('label_',type, sep='')}}])
  for (test in names(results)) {
    result[test] <- results[[test]][,type]
  }
  return(result)
}

ts_results <- .get_results_by_type(results, 'ts')
uci_results <- .get_results_by_type(results, 'uci')
ci_results <- .get_results_by_type(results, 'ci')
lcd_min_results <- .get_results_by_type(results, 'lcd_min')
lcd_CY_results <- .get_results_by_type(results, 'lcd_CY')
CY_results <- .get_results_by_type(results, 'CY')

t1 <- TeX('$(p_{CX} < \\alpha) \\and (p_{XY} < \\alpha) \\and (p_{CY|X} > \\alpha)$')
t2 <- TeX('$(p_{CX} < \\alpha_0) \\and (p_{XY} < \\alpha_0) \\and (p_{CY|X} > \\alpha_0) \\and (p_{CY} < \\alpha)$')
t3 <- TeX('$(p_{CX} < \\alpha) \\and (p_{XY} < \\alpha) \\and (p_{CY|X} > \\min(\\alpha, \\alpha_0))$')
t4 <- TeX('$(p_{CY} < \\alpha) \\and (p_{XY} < \\alpha) \\and (p_{CY|X} > \\min(\\alpha, \\alpha_0))$')
t5 <- TeX('$(p_{CX} < \\alpha) \\and (p_{XY} < \\alpha) \\and (p_{CY|X} > \\\\alpha_0)$')
t6 <- TeX('$(p_{CY} < \\alpha) \\and (p_{XY} < \\alpha) \\and (p_{CY|X} > \\\\alpha_0)$')
grid <- plot_grid(
  pplot_roc(lcd_min_results[,1], lcd_min_results[,-1], t1),
  pplot_roc(lcd_CY_results[,1], lcd_CY_results[,-1], t2),
  pplot_roc_custom(lcd_min_results[,1], ts_results[,-1], uci_results[,-1], ci_results[,-1], t3),
  pplot_roc_custom(lcd_min_results[,1], CY_results[,-1], uci_results[,-1], ci_results[,-1], t4),
  pplot_roc_custom(lcd_min_results[,1], ts_results[,-1], uci_results[,-1], ci_results[,-1], t5, 0),
  pplot_roc_custom(lcd_min_results[,1], CY_results[,-1], uci_results[,-1], ci_results[,-1], t6, 0),
  nrow=1
)
plot(grid)

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
.path <- 'experiments/simulations/output/lcd-measures-roc-tests/'

save.image(file=paste(.path, 'lcd-roc-tests_', timestamp, '.Rdata', sep=''))

.ggsave(paste(.path, timestamp, sep=''), grid, 60, 10)
.ggsave(paste(.path, 'last', sep=''), grid, 60, 10)

