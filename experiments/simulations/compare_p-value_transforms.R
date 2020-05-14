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
  # variance_shift,
  fixed_point
  # mixture
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
      Z <- link_nonlin2 * L
      Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      
      intervene <- rbinom(1, 1, p_link)
      Z <- intervene * do_intervention(interv_options, Z, C) + (1-intervene) * Z
      
      link_nonlin <- link_nonlin1 & link_nonlin2
    }
  }
  
  cond_indep <- as.numeric(cond_indep | !link_nonlin | !intervene)
  
  return(list(C=C, X=X, Z=Z,
              label_ts=as.numeric(intervene),
              label_uci=as.numeric(link_nonlin),
              label_ci=cond_indep))
}

get_results <- function(dataset, test){
  result <- foreach(i=1:length(dataset), .combine=rbind) %dopar% {
    data <- dataset[[i]]
    
    if (test == 'pcor' || test == 'pcor-log') {
      f <- .pcor_wrapper
    } else if (test == 'bayes') {
      f <- .bayes_wrapper
    } else if (test == 'bcor-approx') {
      f <- .bcor_approx_wrapper
    } else {
      f <- .bcor_wrapper
    }
    
    ts <- f(data$C, data$Z)
    uci <- f(data$Z, data$X)
    ci <- f(data$C, data$X, data$Z)
    CX <- f(data$C, data$X)
    
    label <- as.numeric(data$label_uci & data$label_ts & data$label_ci)
    n <- length(data$C)
    
    if (test == 'pcor') {
      return(data.frame(label=label,
                        cx_cond=(ts < 1/sqrt(2*n)) * (uci < 1/sqrt(2*n)) * (ci >= 1/sqrt(2*n)) * (1-CX),
                        ci_cond=(ts < 1/sqrt(2*n)) * (uci < 1/sqrt(2*n)) * ci,
                        min=min(1 - ts, 1 - uci, ci)))
    } else if (test == 'pcor-log') {
      return(data.frame(label=-log(1-label),
                        log_cx_cond=(ts < 1/sqrt(2*n)) * (uci < 1/sqrt(2*n)) * (ci >= 1/sqrt(2*n)) * (-log(CX)),
                        log_ci_cond=(ts < 1/sqrt(2*n)) * (uci < 1/sqrt(2*n)) * (-log(1-ci)),
                        log_min=-log(1-min(1 - ts, 1 - uci, ci))))
    } else if (test == 'bayes' || test == 'bcor' || test == 'bcor-approx') {
      return(data.frame(label=label,
                        cx_cond=(ts < 1/2) * (uci < 1/2) * (ci >= 1/2) * (1-CX),
                        ci_cond=(ts < 1/2) * (uci < 1/2) * ci,
                        min=min(1 - ts, 1 - uci, ci)))
    }
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
pcor_results <- get_results(data, 'pcor')
pcor_log_results <- get_results(data, 'pcor-log')
bcor_results <- get_results(data, 'bcor')
bcor_approx_results <- get_results(data, 'bcor-approx')
bayes_results <- get_results(data, 'bayes')

stopCluster(.cl)


# Process results
##############################################

.pcor_plot <- pplot_roc(pcor_results[,1], pcor_results[,-1], 'pcor')
.pcor_log_plot <- pplot_roc(pcor_log_results[,1], pcor_log_results[,-1], 'pcor-log')
.bcor_plot <- pplot_roc(bcor_results[,1], bcor_results[,-1], 'bcor')
.bcor_approx_plot <- pplot_roc(bcor_approx_results[,1], bcor_approx_results[,-1], 'bcor-approx')
.bayes_plot <- pplot_roc(bayes_results[,1], bayes_results[,-1], 'bayes')

grid <- plot_grid(.pcor_plot, .bcor_approx_plot, .bcor_plot, .bayes_plot, nrow=1)
plot(grid)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save.image(file=paste('experiments/simulations/output/pval-tests/pval_test_',
                      timestamp, ".Rdata", sep=""))

ggsave(
  paste('experiments/simulations/output/pval-tests/pval_test_', timestamp, ".pdf", sep=""),
  plot = grid,
  scale = 1,
  width = 35,
  height = 10,
  units = "cm",
  dpi = 300,
  limitsize = TRUE
)
