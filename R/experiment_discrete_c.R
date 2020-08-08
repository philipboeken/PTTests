experiment_discrete_c <- function (m = 1000, n = -1, graph_probs = c(3/5, 1/5, 1/5), 
                                   dim_C = 4, err_sd = 1/2, p_link = 4/5, 
                                   interv_options = c(mean_shift, variance_shift, 
                                                      fixed_point, mixture),
                                   nonlin_options = c(linear, parabolic, sinusoidal),
                                   seed = 0,
                                   path = 'output/discrete_c/') {
  
  # Setup test
  ##############################################
  set.seed(seed)
  
  n <- if (n < 0) 60 * dim_C else n
  
  get_results <- function(dataset, test){
    `%dopar%` <- foreach::`%dopar%`
    result <- foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
      data <- dataset[[i]]
      
      start_time_CX <- Sys.time()
      CX <- test(data$C, data$X)
      end_time_CX <- Sys.time()
      
      start_time_XY <- Sys.time()
      XY <- test(data$X, data$Y)
      end_time_XY <- Sys.time()
      
      start_time_CY_X <- Sys.time()
      CY_X <- test(data$C, data$Y, data$X)
      end_time_CY_X <- Sys.time()
      
      return(data.frame(
        label_CX = data$label_CX,
        label_XY = data$label_XY,
        label_CY_X = data$label_CY_X,
        label_lcd = data$label_lcd,
        CX = CX, XY = XY, CY_X = CY_X,
        time_CX = end_time_CX - start_time_CX,
        time_XY = end_time_XY - start_time_XY,
        time_CY_X = end_time_CY_X - start_time_CY_X
      ))
    }
    return(result)
  }
  
  .polyatree_wrapper_continuous <- function (X, Y, Z = NULL) {
    if (is.null(Z)) {
      return(polyatree_continuous_independence_test(X, Y)$p_H0)
    }
    
    return(polyatree_continuous_ci_test(X, Y, Z)$p_H0)
  }
  
  # Do test
  ##############################################
  
  doParallel::registerDoParallel()
  
  data <- lapply(1:m, function (i) get_data(graph_probs, n, dim_C, err_sd, p_link, 
                                            interv_options, nonlin_options))
  
  results <- list(
    # ppcor = get_results(data, .ppcor_wrapper),
    rcot = get_results(data, .rcot_wrapper),
    # polyatree_c = get_results(data, .polyatree_wrapper_continuous),
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
