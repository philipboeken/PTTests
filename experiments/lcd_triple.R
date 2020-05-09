source('experiments/test_helpers.R')
library(foreach)
library(doParallel)
library(cowplot)


# Input parameters
##############################################
n <- 300
m <- 500

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
    # lcd <- min((1 - ts), (1 - uci), ci)
    lcd <- (1 - ts) * (1 - uci) * ci
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
  ccit=get_results(data, .ccit_wrapper),
  rcot=get_results(data, .rcot_wrapper)
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

.ts_plot <- pplot_roc(ts_results$label, subset(ts_results, select=-c(label)), 
                      'Two-sample test')
.uci_plot <- pplot_roc(uci_results$label, subset(uci_results, select=-c(label)), 
                       'Continuous independence test')
.ci_plot <- pplot_roc(ci_results$label, subset(ci_results, select=-c(label)), 
                      'Conditional two-sample test')
.lcd_plot <- pplot_roc(lcd_results$label, subset(lcd_results, select=-c(label)), 
                       'LCD ensemble')

grid <- plot_grid(.ts_plot, .uci_plot, .ci_plot, .lcd_plot, nrow=2)
plot(grid)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save.image(file=paste('experiments/output/lcd-roc-tests_', timestamp, ".Rdata", sep=""))

ggsave(
  paste('experiments/output/lcd-roc-tests_', timestamp, ".pdf", sep=""),
  plot = grid,
  scale = 1,
  width = 20,
  height = 20,
  units = "cm",
  dpi = 300,
  limitsize = TRUE
)
