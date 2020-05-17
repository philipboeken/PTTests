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
n <- 400
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
    ts <- test(data$C, data$Z)
    uci <- test(data$Z, data$X)
    ci <- test(data$C, data$X, data$Z)
    lcd <- min((1 - ts), (1 - uci), ci)
    label_lcd <- as.numeric(data$label_uci & data$label_ts & data$label_ci)
    return(data.frame(
      label_ts=data$label_ts,
      label_uci=data$label_uci,
      label_ci=data$label_ci,
      label_lcd=label_lcd,
      ts=1-ts,
      uci=1-uci,
      ci=ci,
      lcd=lcd
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
  pcor=get_results(data, .pcor_wrapper),
  pcor_bayes=get_results(data, .bayes_transform(.pcor_wrapper)),
  bayes=get_results(data, .bayes_wrapper),
  bcor_wg=get_results(data, .bcor_wg_wrapper),
  # bcor_approx=get_results(data, .bcor_approx_wrapper),
  bcor_ly=get_results(data, .bcor_ly_wrapper),
  gcm_bayes=get_results(data, .bayes_transform(.gcm_wrapper)),
  # gcm=get_results(data, .gcm_wrapper),
  # ccit=get_results(data, .ccit_wrapper),
  rcot_bayes=get_results(data, .bayes_transform(.rcot_wrapper))
  # rcot=get_results(data, .rcot_wrapper)
)

stopCluster(.cl)


# Process results
##############################################

.get_results_by_type <- function (results, type) {
  result <- data.frame(label=results$bayes[,{{paste('label_',type, sep="")}}])
  for (test in names(results)) {
    result[test] <- results[[test]][,type]
  }
  return(result)
}

ts_results <- .get_results_by_type(results, 'ts')
uci_results <- .get_results_by_type(results, 'uci')
ci_results <- .get_results_by_type(results, 'ci')
lcd_results <- .get_results_by_type(results, 'lcd')

grid <- plot_grid(
  pplot_roc(ts_results[,1], ts_results[,-1], 'Two-sample test'), 
  pplot_roc(uci_results[,1], uci_results[,-1], 'Continuous independence test'), 
  pplot_roc(ci_results[,1], ci_results[,-1], 'Conditional two-sample test'), 
  pplot_roc(lcd_results[,1], lcd_results[,-1], 'LCD ensemble'), 
  nrow=2
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
.path <- 'experiments/simulations/output/lcd-roc-tests/'

save.image(file=paste(.path, 'lcd-roc-tests_', timestamp, ".Rdata", sep=""))

.ggsave(paste(.path, 'lcd-roc-tests_', timestamp, sep=""), grid, 20, 20)
.ggsave(paste(.path, 'lcd-roc-tests_last', sep=""), grid, 20, 20)

# plot(grid)
