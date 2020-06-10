experiment_lcd_compare_tests <- function (n = 400, m = 2000, err_sd = 0.5, 
                                          p_ci = 0.6, p_link = 0.8, p_two_sample = 0.5 ,
                                          nonlin_options = c(.linear, .parabolic, .sinusoidal),
                                          interv_options = c(.mean_shift, .variance_shift, 
                                                             .fixed_point, .mixture),
                                          seed = 0,
                                          path = 'output/lcd_compare_tests/') {
  
  # Setup test
  ##############################################
  set.seed(seed)
  
  get_data <- function() {
    C <- rbinom(n, 1, p_two_sample)
    
    cond_indep <- rbinom(1, 1, p_ci)
    if (cond_indep) { # C -> X -> Y
      intervene <- rbinom(1, 1, p_link)
      X <- intervene * .do_intervention(interv_options, rnorm(n), C) + (1-intervene) * rnorm(n)
      
      link_nonlin <- rbinom(1, 1, p_link)
      Y <- link_nonlin * .nonlin(nonlin_options, X)
      Y <- Y + err_sd * rnorm(n, 0, ifelse(sd(Y) > 0, sd(Y), 1/err_sd))
    } else {
      if (runif(1) <= 0.5) { # C -> X <- Y
        Y <- rnorm(n)
        
        link_nonlin <- rbinom(1, 1, p_link)
        X <- link_nonlin * .nonlin(nonlin_options, Y)
        X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1/err_sd))
        
        intervene <- rbinom(1, 1, p_link)
        X <- intervene * .do_intervention(interv_options, X, C) + (1-intervene) * X
      } else { # C -> X <- L -> Y
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
    `%dopar%` <- foreach::`%dopar%`
    result <- foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
      data <- dataset[[i]]
      
      start_time_ts <- Sys.time()
      ts <- test(data$C, data$X)
      end_time_ts <- Sys.time()
      
      start_time_uci <- Sys.time()
      uci <- test(data$X, data$Y)
      end_time_uci <- Sys.time()
      
      start_time_ci <- Sys.time()
      ci <- test(data$C, data$Y, data$X)
      end_time_ci <- Sys.time()
      
      return(data.frame(
        label_ts = data$label_ts,
        label_uci = data$label_uci,
        label_ci = data$label_ci,
        label_lcd = data$label_lcd,
        ts = ts,
        uci = uci,
        ci = ci,
        time_ts = end_time_ts - start_time_ts,
        time_uci = end_time_uci - start_time_uci,
        time_ci = end_time_ci - start_time_ci
      ))
    }
    return(result)
  }
  
  # Do test
  ##############################################
  cores <- parallel::detectCores()
  cl <- parallel::makeForkCluster(cores[1]-1)
  doParallel::registerDoParallel(cl)
  
  data <- lapply(1:m, function (i) get_data())
  results <- list(
    ppcor = get_results(data, .ppcor_wrapper),
    spcor = get_results(data, .spcor_wrapper),
    ppcor_b = get_results(data, .ppcor_b_wrapper),
    gcm = get_results(data, .gcm_wrapper),
    rcot = get_results(data, .rcot_wrapper),
    # ccit = get_results(data, .ccit_wrapper),
    polyatree = get_results(data, .polyatree_wrapper)
  )
  
  parallel::stopCluster(cl)
  
  
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
  
  `%do%` <- foreach::`%do%`
  times <- foreach::foreach(test = names(results), .combine = rbind) %do% {
    rbind(c(test, '1_ts', sum(results[[test]][,'time_ts'])),
          c(test, '2_uci', sum(results[[test]][,'time_uci'])),
          c(test, '3_ci', sum(results[[test]][,'time_ci'])))
  }
  times <- data.frame(ensemble = times[,1], test = times[,2], time = as.double(times[,3]))
  
  ts_results <- .get_results_by_type(results, 'ts')
  uci_results <- .get_results_by_type(results, 'uci')
  ci_results <- .get_results_by_type(results, 'ci')
  
  labels_lcd <- factor(results$polyatree[,'label_lcd'], ordered = TRUE, levels = c(1,0))
  
  plot_ts <- .plot_roc(ts_results[,1], ts_results[,-1], freq_default = 0.05)
  plot_uci <- .plot_roc(uci_results[,1], uci_results[,-1], freq_default = 0.05)
  plot_ci <- .plot_roc(ci_results[,1], ci_results[,-1], freq_default = 0.05)
  plot_lcd <- .plot_roc_custom(labels_lcd, ts_results[,-1], uci_results[,-1], ci_results[,-1])
  
  grid <- cowplot::plot_grid(plot_ts, plot_uci, plot_ci, plot_lcd, nrow = 1)
  
  plot_runtimes <- .plot_times(times)
  
  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
  
  save(results, file = sprintf("%s%s.Rdata", path, timestamp))
  .ggsave(paste(path, timestamp, sep = ''), grid, 40, 10)
  .ggsave(paste(path, 'last', sep = ''), grid, 40, 10)
  .ggsave(paste(path, 'two-sample', sep = ''), plot_ts, 10, 10)
  .ggsave(paste(path, 'marginal-independence', sep = ''), plot_uci, 10, 10)
  .ggsave(paste(path, 'conditional-independence', sep = ''), plot_ci, 10, 10)
  .ggsave(paste(path, 'lcd', sep = ''), plot_lcd, 10, 10)
  .ggsave(paste(path, 'times', sep = ''), plot_runtimes, 20, 8)
}
