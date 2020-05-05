source('init.R', chdir=TRUE)
source('experiments/test_wrappers.R')
library(foreach)
library(doParallel)
library(cowplot)

# Define mappings
##############################################
mean_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (base + theta)
}

variance_shift <- function(base, C) {
  theta <- sample(1:3, 1)
  (1-C) * base + C * (1+theta) * base
}

mixture <- function(base, C) {
  theta <- sample(1:3, 1)
  idx <- sample(c(-1,1), length(C), replace = TRUE)
  (1-C) * base + C * (base + idx*theta)
}

do_intervention <- function (base, C) {
  mapping <- sample(c(
    mean_shift,
    variance_shift,
    mixture
  ), 1)[[1]]
  return(mapping(base, C))
}

linear <- function(X) {
  2*X / 3
}

parabolic <- function(X) {
  2*X^2 / 3
}

sinusoidal <- function(X) {
  2*sin(3*X)
}

partial <- function(X) {
  b <- rbinom(length(X), 1, 0.1)
  b*X + (1-b)*rnorm(length(X))
}

nonlin <- function(X) {
  mapping <- sample(c(
    linear,
    parabolic,
    sinusoidal,
    partial
  ), 1)[[1]]
  return(mapping(X))
}


# Setup test
##############################################
get_data <- function(n) {
  p <- runif(1, 0.45, 0.65)
  C <- rbinom(n, 1, p)
  
  cond_indep <- sample(c(0, 1), 1)
  if (cond_indep) { # C -> Z -> X
    intervene <- sample(c(0,1), 1)
    Z <- intervene * do_intervention(rnorm(n), C) + (1-intervene) * rnorm(n)
    
    link_nonlin <- sample(c(0,1), 1)
    errs <- rnorm(n, 0, runif(1, 1, 2.5))
    X <- link_nonlin * nonlin(Z) + errs
  } else {
    if (runif(1) <= 0) { # C -> Z <- X
      X <- rnorm(n)
      
      link_nonlin <- sample(c(0,1), 1)
      errs <- rnorm(n, 0, runif(1, 1, 1.5))
      Z <- link_nonlin * nonlin(X) + errs
      
      intervene <- sample(c(0,1), 1)
      Z <- intervene * do_intervention(Z, C) + (1-intervene) * Z
    } else { # C -> Z <- L -> X
      L <- rnorm(n)
      
      link_nonlin1 <- sample(c(0,1), 1)
      errs <- rnorm(n, 0, runif(1, 1, 2.5))
      X <- link_nonlin1 * nonlin(L) + errs
      
      link_nonlin2 <- sample(c(0,1), 1)
      errs <- rnorm(n, 0, runif(1, 1, 2.5))
      Z <- link_nonlin2 * nonlin(L) + errs
      
      intervene <- sample(c(0,1), 1)
      Z <- intervene * do_intervention(Z, C) + (1-intervene) * Z
      
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
  result <- foreach(i=1:length(dataset), .combine=rbind) %do% {
    data <- dataset[[i]]
    ci <- test(data$X, data$C, data$Z)
    uci <- test(data$X, data$Z)
    ts <- test(data$Z, data$C)
    lcd <- min((1 - ts), (1 - uci), ci)
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
cores <- detectCores()
cl <- makeForkCluster(cores[1]-1)
registerDoParallel(cl)
n <- 400
m <- 400
data <- lapply(1:m, function (i) get_data(n))
results <- list(
  # pcor=get_results(data, pcor_wrapper)
  bayes=get_results(data, bayes_ci_wrapper),
  # gcm=get_results(data, gcm_wrapper)
  rcot=get_results(data, rcot_wrapper)
)

stopCluster(cl)


# Process results
##############################################
uci_results <- data.frame(
  # pcor=results$pcor[,'uci'],
  bayes=results$bayes[,'uci'],
  # gcm=results$gcm[,'uci'],
  rcot=results$rcot[,'uci']
)
uci_label <- results$bayes[,'label_uci']

ci_results <- data.frame(
  # pcor=results$pcor[,'ci'],
  bayes=results$bayes[,'ci'],
  # gcm=results$gcm[,'ci'],
  rcot=results$rcot[,'ci']
)
ci_label <- results$bayes[,'label_ci']

ts_results <- data.frame(
  # pcor=results$pcor[,'ts'],
  bayes=results$bayes[,'ts'],
  # gcm=results$gcm[,'ts'],
  rcot=results$rcot[,'ts']
)
ts_label <- results$bayes[,'label_ts']

lcd_results <- data.frame(
  # pcor=results$pcor[,'lcd'],
  bayes=results$bayes[,'lcd'],
  # gcm=results$gcm[,'lcd'],
  rcot=results$rcot[,'lcd']
)
lcd_label <- results$bayes[,'label_lcd']

uci_plot <- pplot_roc(uci_label, uci_results, 'X_||_Z')
ci_plot <- pplot_roc(ci_label, ci_results, 'X_||_C|Z')
ts_plot <- pplot_roc(ts_label, ts_results, 'Z_||_C')
lcd_plot <- pplot_roc(lcd_label, lcd_results, 'LCD: Z -> X')

grid <- plot_grid(ts_plot, uci_plot, ci_plot, lcd_plot, nrow=2)
plot(grid)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save.image(file=paste('experiments/output/lcd-roc-tests_', timestamp, ".Rdata", sep=""))

ggsave(
  paste('experiments/output/lcd-roc-tests_', timestamp, ".pdf", sep=""),
  plot = grid,
  scale = 1,
  width = 30,
  height = 15,
  # height = 7.5,
  units = "cm",
  dpi = 300,
  limitsize = TRUE
)


