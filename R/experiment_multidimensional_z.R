experiment_multidimensional_z <- function (m = 1000, n = -1, dim_C = 2, err_sd = 1/2, p_link = 4/5, 
                                           interv_options = c(mean_shift, variance_shift, 
                                                              fixed_point, mixture),
                                           nonlin_options = c(linear, parabolic, sinusoidal),
                                           seed = 0,
                                           path = 'output/multidimensional_z/') {
  
  # Setup test
  ##############################################
  set.seed(seed)
  
  n <- if (n < 0) 40 * dim_C else n
  
  get_data <- function(n, dim_C, dim_X, err_sd, p_link) {
    # X_0 - X_1 - ... - X_{dim_X} - X_{dim_X+1}
    
    link_nonlin <- sample(c(-1, 0, 1), dim_X+1, replace = TRUE, prob = c(p_link/2, 1-p_link, p_link/2))
    X <- matrix(0, ncol=dim_X+2, nrow=n)
    X[, which(diff(c(-1, link_nonlin, 1)) %in% c(1, 2))] <- rnorm(n)
    updated <- which(diff(c(-1, link_nonlin, 1)) %in% c(1, 2))
    
    to_update <- setdiff(1:(dim_X+2), updated)
    
    while(length(to_update) > 0) {
      child <- if (length(to_update > 1)) sample(to_update, 1) else to_update
      all_parents <- setdiff(c(ifelse(link_nonlin[child-1] == 1, child-1, 0), 
                               ifelse(link_nonlin[child] == -1, child+1, 0)), 0)
      updated_parents <- setdiff(all_parents, to_update)
      for (parent in updated_parents) {
        Z <- .nonlin(nonlin_options, X[, parent])
        X[, child] <- X[, child] + Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
      }
      if (length(all_parents) == length(updated_parents) && all(all_parents == updated_parents)) {
        to_update <- setdiff(to_update, child)
      }
    }
    
    colliders <- 1+which(diff(link_nonlin) == -2)
    non_colliders <- setdiff(2:(2+dim_X), colliders)
    
    return(list(X = X, label = length(non_colliders) > 0))
  }
  
  get_results <- function(dataset, test, dim_X){
    `%dopar%` <- foreach::`%dopar%`
    result <- foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
      data <- dataset[[i]]
      X <- data$X[,-c(1, dim_X+2)]
      
      start_time <- Sys.time()
      result <- test(data$C, data$Y, data$X)
      end_time <- Sys.time()
      
      return(data.frame(
        label = data$label,
        result = result,
        time = end_time - start_time
      ))
    }
    return(result)
  }
  
  # Do test
  ##############################################
  
  doParallel::registerDoParallel()
  
  data <- lapply(1:m, function (i) get_data(n, dim_C, err_sd, p_link))
  
  results <- list(
    # ppcor = get_results(data, .ppcor_wrapper),
    rcot = get_results(data, .rcot_wrapper),
    polyatree = get_results(data, .polyatree_wrapper)
  )
  
  
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
    rbind(c(test, '1_CX', sum(results[[test]][,'time_CX'])),
          c(test, '2_XY', sum(results[[test]][,'time_XY'])),
          c(test, '3_CY_X', sum(results[[test]][,'time_CY_X'])))
  }
  times <- data.frame(ensemble = times[,1], test = times[,2], time = as.double(times[,3]))
  
  CX_results <- .get_results_by_type(results, 'CX')
  XY_results <- .get_results_by_type(results, 'XY')
  CY_X_results <- .get_results_by_type(results, 'CY_X')
  
  labels_lcd <- factor(results$polyatree[,'label_lcd'], ordered = TRUE, levels = c(1,0))
  
  plot_CX <- .plot_roc(CX_results[,1], CX_results[,-1], freq_default = 0.05)
  plot_XY <- .plot_roc(XY_results[,1], XY_results[,-1], freq_default = 0.05)
  plot_CY_X <- .plot_roc(CY_X_results[,1], CY_X_results[,-1], freq_default = 0.05)
  plot_lcd <- .plot_roc_custom(labels_lcd, CX_results[,-1], XY_results[,-1], CY_X_results[,-1])
  
  grid <- cowplot::plot_grid(plot_CX, plot_XY, plot_CY_X, plot_lcd, nrow = 1)
  
  plot_runtimes <- .plot_times(times)
  
  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
  
  save(results, file = sprintf("%s%s.Rdata", path, timestamp))
  .ggsave(paste(path, timestamp, sep = ''), grid, 40, 10)
  .ggsave(paste(path, 'last', sep = ''), grid, 40, 10)
  .ggsave(paste(path, 'two-sample', sep = ''), plot_CX, 10, 10)
  .ggsave(paste(path, 'marginal-independence', sep = ''), plot_XY, 10, 10)
  .ggsave(paste(path, 'conditional-independence', sep = ''), plot_CY_X, 10, 10)
  .ggsave(paste(path, 'lcd', sep = ''), plot_lcd, 10, 10)
  .ggsave(paste(path, 'times', sep = ''), plot_runtimes, 20, 8)
}
