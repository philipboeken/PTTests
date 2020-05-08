source('experiments/test_helpers.R')
library(foreach)
library(doParallel)
library(cowplot)


# Input parameters
##############################################
n <- 300
m <- 400

p_link <- 0.9

err_sd <- 0.1
nonlin_options <- c(
  linear,
  parabolic,
  sinusoidal,
  partial
)

p_two_sample <- 0.5
interv_options <- c(
  mean_shift,
  variance_shift,
  mixture
  # tails
)

p_ci <- 0.8


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
  
  return(list(C=C, Z=Z, X=X,
              label_ci=cond_indep,
              label_ts=as.numeric(!intervene),
              label_uci=as.numeric(!link_nonlin)))
}

get_results <- function(dataset, test){
  result <- foreach(i=1:length(dataset), .combine=rbind) %dopar% {
    data <- dataset[[i]]
    ci <- test(data$X, data$C, data$Z)
    uci <- test(data$X, data$Z)
    ts <- test(data$Z, data$C)
    lcd <- min((1 - ts), (1 - uci), ci)
    # lcd <- (1 - ts) * (1 - uci) * ci
    label_lcd <- as.numeric(!data$label_uci & !data$label_ts & data$label_ci)
    return(data.frame(label_ci=data$label_ci,
                      label_uci=data$label_uci,
                      label_ts=data$label_ts,
                      label_lcd=label_lcd,
                      ci=ci,
                      uci=uci,
                      ts=ts,
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
  bayes=get_results(data, .bayes_wrapper),
  gcm=get_results(data, .gcm_wrapper),
  rcot=get_results(data, .rcot_wrapper),
  ccit=get_results(data, .ccit_wrapper)
)

stopCluster(.cl)


# Process results
##############################################
uci_results <- data.frame(
  pcor=results$pcor[,'uci'],
  bayes=results$bayes[,'uci'],
  gcm=results$gcm[,'uci'],
  rcot=results$rcot[,'uci'],
  ccit=results$ccit[,'uci']
)
uci_label <- results$bayes[,'label_uci']

ci_results <- data.frame(
  pcor=results$pcor[,'ci'],
  bayes=results$bayes[,'ci'],
  gcm=results$gcm[,'ci'],
  rcot=results$rcot[,'ci'],
  ccit=results$ccit[,'ci']
)
ci_label <- results$bayes[,'label_ci']

ts_results <- data.frame(
  pcor=results$pcor[,'ts'],
  bayes=results$bayes[,'ts'],
  gcm=results$gcm[,'ts'],
  rcot=results$rcot[,'ts'],
  ccit=results$ccit[,'ts']
)
ts_label <- results$bayes[,'label_ts']

lcd_results <- data.frame(
  pcor=results$pcor[,'lcd'],
  bayes=results$bayes[,'lcd'],
  gcm=results$gcm[,'lcd'],
  rcot=results$rcot[,'lcd'],
  ccit=results$ccit[,'lcd']
)
lcd_label <- results$bayes[,'label_lcd']

ts_plot <- pplot_roc(ts_label, ts_results, 'Two-sample test')
uci_plot <- pplot_roc(uci_label, uci_results, 'Continuous independence test')
ci_plot <- pplot_roc(ci_label, ci_results, 'Conditional two-sample test')
lcd_plot <- pplot_roc(lcd_label, lcd_results, 'LCD ensemble')

grid <- plot_grid(ts_plot, uci_plot, ci_plot, lcd_plot, nrow=2)
plot(grid)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save.image(file=paste('experiments/output/lcd-roc-tests_', timestamp, ".Rdata", sep=""))

ggsave(
  paste('experiments/output/lcd-roc-tests_', timestamp, ".pdf", sep=""),
  plot = grid,
  scale = 1,
  width = 20,
  height = 20,
  # height = 7.5,
  units = "cm",
  dpi = 300,
  limitsize = TRUE
)
