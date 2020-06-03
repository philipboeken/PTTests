# Imports
source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('helpers.R')
library(foreach)
library(doParallel)
library(cowplot)


# Input parameters
##############################################
set.seed(0)

n <- 400
m <- 2000

err_sd <- 0.5

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
get_data <- function() {
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
              label_ts=1-as.numeric(intervene),
              label_uci=1-as.numeric(link_nonlin),
              label_ci=cond_indep,
              label_lcd=lcd))
}

get_results <- function(dataset, test){
  result <- foreach(i=1:length(dataset), .combine=rbind) %dopar% {
    data <- dataset[[i]]
    
    start_time_ts <- Sys.time()
    ts <- test(data$C, data$Z)
    end_time_ts <- Sys.time()
    
    start_time_uci <- Sys.time()
    uci <- test(data$Z, data$X)
    end_time_uci <- Sys.time()
    
    start_time_ci <- Sys.time()
    ci <- test(data$C, data$X, data$Z)
    end_time_ci <- Sys.time()
    
    return(data.frame(
      label_ts=data$label_ts,
      label_uci=data$label_uci,
      label_ci=data$label_ci,
      label_lcd=data$label_lcd,
      ts=ts,
      uci=uci,
      ci=ci,
      time_ts=end_time_ts - start_time_ts,
      time_uci=end_time_uci - start_time_uci,
      time_ci=end_time_ci - start_time_ci
    ))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)
data <- lapply(1:m, function (i) get_data())
results <- list(
  ppcor=get_results(data, .pcor_wrapper),
  spcor=get_results(data, .prcor_wrapper),
  ppcor_b=get_results(data, .bcor_wrapper),
  gcm=get_results(data, .gcm_wrapper),
  rcot=get_results(data, .rcot_wrapper),
  ccit=get_results(data, .ccit_wrapper),
  polyatree=get_results(data, .bayes_wrapper)
)

stopCluster(.cl)


# Process results
##############################################

.get_results_by_type <- function (results, type) {
  labels <- results$polyatree[,{{paste('label_',type, sep='')}}]
  labels <- factor(labels, ordered = TRUE, levels = c(1, 0))
  result <- data.frame(label=labels)
  for (test in names(results)) {
    result[test] <- results[[test]][,type]
  }
  return(result)
}

times <- foreach(test=names(results), .combine=rbind) %do% {
  rbind(c(test, '1_ts', sum(results[[test]][,'time_ts'])),
        c(test, '2_uci', sum(results[[test]][,'time_uci'])),
        c(test, '3_ci', sum(results[[test]][,'time_ci'])))
}
times <- data.frame(ensemble=times[,1], test=times[,2], time=as.double(times[,3]))

ts_results <- .get_results_by_type(results, 'ts')
uci_results <- .get_results_by_type(results, 'uci')
ci_results <- .get_results_by_type(results, 'ci')

labels_lcd <- factor(results$polyatree[,'label_lcd'], ordered = TRUE, levels = c(1,0))

.plot_ts <- plot_roc(ts_results[,1], ts_results[,-1], freq_default=0.05)
.plot_uci <- plot_roc(uci_results[,1], uci_results[,-1], freq_default=0.05)
.plot_ci <- plot_roc(ci_results[,1], ci_results[,-1], freq_default=0.05)
.plot_lcd <- plot_roc_custom(labels_lcd, ts_results[,-1], uci_results[,-1], ci_results[,-1])

grid <- plot_grid(.plot_ts, .plot_uci, .plot_ci, .plot_lcd, nrow=1)

.plot_time <- plot_times(times)
  
timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
.path <- 'experiments/simulations/output/lcd-roc-tests/'

save.image(file=paste(.path, timestamp, '.Rdata', sep=''))

.ggsave(paste(.path, timestamp, sep=''), grid, 40, 10)
.ggsave(paste(.path, 'last', sep=''), grid, 40, 10)
.ggsave(paste(.path, 'two-sample', sep=''), .plot_ts, 10, 10)
.ggsave(paste(.path, 'marginal-independence', sep=''), .plot_uci, 10, 10)
.ggsave(paste(.path, 'conditional-independence', sep=''), .plot_ci, 10, 10)
.ggsave(paste(.path, 'lcd', sep=''), .plot_lcd, 10, 10)
.ggsave(paste(.path, 'times', sep=''), .plot_time, 20, 8)
