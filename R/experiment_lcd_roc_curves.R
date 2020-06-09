experiment_lcd_roc_curves <- function (n = 400, m = 2000, err_sd = 0.5, 
                                       p_ci = 0.6, p_link = 0.8, p_two_sample = 0.5 ,
                                       .nonlin_options = c(.linear, .parabolic, .sinusoidal),
                                       interv_options = c(.mean_shift, .variance_shift, 
                                                          .fixed_point, .mixture),
                                       seed = 0,
                                       path = 'output/lcd_roc_curves/') {
  library(doParallel)
  
  set.seed(seed)
  
  # Setup test
  ##############################################
  get_data <- function() {
    C <- rbinom(n, 1, p_two_sample)
    
    cond_indep <- rbinom(1, 1, p_ci)
    if (cond_indep) { # C -> X -> Y
      intervene <- rbinom(1, 1, p_link)
      X <- intervene * .do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
      
      link_nonlin <- rbinom(1, 1, p_link)
      Y <- link_nonlin * .nonlin(.nonlin_options, X)
      Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
    } else {
      if (runif(1) <= 0.5) { # C -> X <- Y
        Y <- rnorm(n)
        
        link_nonlin <- rbinom(1, 1, p_link)
        X <- link_nonlin * .nonlin(.nonlin_options, Y)
        X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
        
        intervene <- rbinom(1, 1, p_link)
        X <- intervene * .do_intervention(interv_options, X, C) + (1-intervene) * X
      } else { # C -> X <- L -> Y
        L <- rnorm(n)
        
        link_nonlin1 <- rbinom(1, 1, p_link)
        Y <- link_nonlin1 * .nonlin(.nonlin_options, L)
        Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
        
        link_nonlin2 <- rbinom(1, 1, p_link)
        X <- link_nonlin2 * .nonlin(.nonlin_options, L)
        X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
        
        intervene <- rbinom(1, 1, p_link)
        X <- intervene * .do_intervention(interv_options, X, C) + (1-intervene) * X
        
        link_nonlin <- link_nonlin1 & link_nonlin2
      }
    }
    
    cond_indep <- as.numeric(cond_indep | !link_nonlin | !intervene)
    lcd <- as.numeric(intervene & link_nonlin & cond_indep)
    
    return(list(C = C, X = X, Y = Y,
                label_ts = 1-as.numeric(intervene),
                label_uci = 1-as.numeric(link_nonlin),
                label_ci = cond_indep,
                label_lcd = lcd))
  }
  
  get_results <- function(dataset, test){
    result <- foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
      data <- dataset[[i]]
      
      ts <- test(data$C, data$X)
      uci <- test(data$X, data$Y)
      ci <- test(data$C, data$Y, data$X)
      
      return(data.frame(
        label_ts = data$label_ts,
        label_uci = data$label_uci,
        label_ci = data$label_ci,
        label_lcd = data$label_lcd,
        ts = ts,
        uci = uci,
        ci = ci
      ))
    }
    return(result)
  }
  
  
  # Do test
  ##############################################
  cores <- detectCores()
  cl <- makeForkCluster(cores[1]-1)
  registerDoParallel(cl)
  data <- lapply(1:m, function (i) get_data())
  results <- list(
    ppcor = get_results(data, .ppcor_wrapper),
    ppcor_b = get_results(data, .ppcor_b_wrapper),
    polyatree = get_results(data, .polyatree_wrapper)
  )
  
  stopCluster(cl)
  
  
  # Process results
  ##############################################
  
  .get_results_by_type <- function (results, type) {
    labels <- results$polyatree[,{{paste('label_',type, sep = '')}}]
    labels <- factor(labels, ordered = TRUE, levels = c(1, 0))
    result <- data.frame(label = labels)
    for (test in names(results)) {
      result[test] <- results[[test]][,type]
    }
    return(result)
  }
  
  ts_results <- .get_results_by_type(results, 'ts')
  uci_results <- .get_results_by_type(results, 'uci')
  ci_results <- .get_results_by_type(results, 'ci')
  
  labels_lcd <- factor(results$polyatree[,'label_lcd'], ordered = TRUE, levels = c(1,0))
  
  t0 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > \\min(\\alpha_0, 1-\\alpha))$')
  t1 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > \\alpha)$')
  t2 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > 1-\\alpha)$')
  .plot0 <- .plot_roc_custom(labels_lcd, ts_results[,-1], uci_results[,-1], ci_results[,-1], t0, option = 0, plot_point = FALSE)
  .plot1 <- .plot_roc_custom(labels_lcd, ts_results[,-1], uci_results[,-1], ci_results[,-1], t1, option = 1, plot_point = FALSE)
  .plot2 <- .plot_roc_custom(labels_lcd, ts_results[,-1], uci_results[,-1], ci_results[,-1], t2, option = 2, plot_point = FALSE)
  
  grid <- cowplot::plot_grid(.plot0, .plot1, .plot2, nrow = 1)
  
  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
  
  save(results, file = sprintf("%s%s.Rdata", path, timestamp))
  .ggsave(paste(path, timestamp, sep = ''), grid, 30, 10)
  .ggsave(paste(path, 'last', sep = ''), grid, 30, 10)
  .ggsave(paste(path, 'plot0', sep = ''), .plot0, 10, 10)
  .ggsave(paste(path, 'plot1', sep = ''), .plot1, 10, 10)
  .ggsave(paste(path, 'plot2', sep = ''), .plot2, 10, 10)
}
