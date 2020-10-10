experiment_lcd_compare_tests_roc <- function(m = 2000, n = 400, graph_probs = c(3 / 5, 1 / 5, 1 / 5),
                                         dim_C = 2, err_sd = 1 / 2, p_link = 4 / 5,
                                         interv_options = c(mean_shift, variance_shift, mixture),
                                         nonlin_options = c(linear, parabolic, sinusoidal),
                                         simulation = NULL,
                                         seed = 0,
                                         path = 'output/lcd_compare_tests/',
                                         save_figures = TRUE,
                                         pt_continuous = FALSE) {

  # Setup test
  ##############################################

  set.seed(seed)

  if (!is.null(simulation) && simulation == 'paper') {
    interv_options <- c(mean_shift_paper, variance_shift_paper, fixed_point_paper, mixture_paper)
    nonlin_options <- c(linear_paper, parabolic_paper, sinusoidal_paper)
  } else if (!is.null(simulation) && simulation == 'thesis') {
    interv_options <- c(mean_shift, variance_shift, mixture)
    nonlin_options <- c(linear, parabolic, sinusoidal)
  }

  get_results <- function(dataset, test) {
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

  # Do test
  ##############################################

  doParallel::registerDoParallel()

  if (!is.null(simulation) && simulation == 'paper') {
    data <- lapply(1:m, function(i) get_data_paper(graph_probs, n, dim_C, err_sd, p_link,
                                                   interv_options, nonlin_options))
  } else {
    data <- lapply(1:m, function(i) get_data(graph_probs, n, dim_C, err_sd, p_link,
                                             interv_options, nonlin_options))
  }

  if (pt_continuous) {
    results <- list(
      ppcor = get_results(data, .ppcor_wrapper),
      spcor = get_results(data, .spcor_wrapper),
      gcm = get_results(data, .gcm_wrapper),
      ccit = get_results(data, .ccit_wrapper),
      rcot = get_results(data, .rcot_wrapper),
      polyatree_c = get_results(data, .pt_wrapper_continuous),
      polyatree = get_results(data, .pt_wrapper))
  } else {
    results <- list(
      ppcor = get_results(data, .ppcor_wrapper),
      spcor = get_results(data, .spcor_wrapper),
      gcm = get_results(data, .gcm_wrapper),
      ccit = get_results(data, .ccit_wrapper),
      rcot = get_results(data, .rcot_wrapper),
      polyatree = get_results(data, .pt_wrapper))
  }


  # Process results
  ##############################################

  .get_results_by_type <- function(results, type) {
    labels <- results$polyatree[, {{ paste('label_', type, sep = '') }}]
    labels <- factor(labels, ordered = TRUE, levels = c(1, 0))
    result <- data.frame(label = labels)
    for (test in names(results)) {
      result[test] <- results[[test]][, type]
    }
    return(result)
  }

  `%do%` <- foreach::`%do%`
  times <- foreach::foreach(test = names(results), .combine = rbind) %do% {
    rbind(c(test, '1_CX', sum(results[[test]][, 'time_CX'])),
          c(test, '2_XY', sum(results[[test]][, 'time_XY'])),
          c(test, '3_CY_X', sum(results[[test]][, 'time_CY_X'])))
  }
  times <- data.frame(ensemble = times[, 1], test = times[, 2], time = as.double(times[, 3]))

  CX_results <- .get_results_by_type(results, 'CX')
  XY_results <- .get_results_by_type(results, 'XY')
  CY_X_results <- .get_results_by_type(results, 'CY_X')

  labels_lcd <- factor(results$polyatree[, 'label_lcd'], ordered = TRUE, levels = c(1, 0))

  if (save_figures) {
    plot_CX <- .plot_roc(CX_results[, 1], CX_results[, -1], freq_default = 0.05, legend_pos = "none")
    plot_XY <- .plot_roc(XY_results[, 1], XY_results[, -1], freq_default = 0.05, legend_pos = "none")
    plot_CY_X <- .plot_roc(CY_X_results[, 1], CY_X_results[, -1], freq_default = 0.05, legend_pos = "none")
    plot_lcd <- .plot_roc_lcd(labels_lcd, CX_results[, -1], XY_results[, -1], CY_X_results[, -1], legend_pos = "none")

    grid <- cowplot::plot_grid(plot_CX, plot_XY, plot_CY_X, plot_lcd, nrow = 1)

    plot_runtimes <- .plot_roc_times(times)

    timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

    save(results, file = sprintf("%s%s.Rdata", path, timestamp))
    .ggsave(paste(path, timestamp, sep = ''), grid, 40, 10)
    .ggsave(paste(path, 'last', sep = ''), grid, 40, 10)
    .ggsave(paste(path, 'two-sample', timestamp, sep = ''), plot_CX, 8, 8)
    .ggsave(paste(path, 'marginal-independence', timestamp, sep = ''), plot_XY, 8, 8)
    .ggsave(paste(path, 'conditional-independence', timestamp, sep = ''), plot_CY_X, 8, 8)
    .ggsave(paste(path, 'lcd', timestamp, sep = ''), plot_lcd, 8, 8)
    .ggsave(paste(path, 'times', timestamp, sep = ''), plot_runtimes, 15, 6)
  }

  get_test_aucs <- function(labels, predictions) {
    auc_data <- c()
    predictions <- as.matrix(predictions)
    for (i in 1:ncol(predictions)) {
      name <- colnames(predictions)[i]
      roc <- .get_roc(labels, - predictions[, i])
      auc_data[[name]] <- roc$auc
    }
    return(auc_data)
  }

  get_lcd_aucs <- function(labels, CX, XY, CY_X) {
    auc_data <- c()
    for (i in 1:ncol(CX)) {
      name <- colnames(CX)[i]
      bayes <- (name == 'polyatree' || name == 'ppcor_b' || name == 'polyatree_c')
      roc <- .get_lcd_roc(labels, CX[, i], XY[, i], CY_X[, i], bayes)
      auc_data[[name]] <- roc$auc
    }
    return(auc_data)
  }

  times_trans <- sapply(unique(times[, 'ensemble']),
                        function(test) sum(times[times[, 'ensemble'] == test, 'time']),
                        simplify = FALSE, USE.NAMES = TRUE)

  return(list(
    CX = get_test_aucs(CX_results[, 1], CX_results[, -1]),
    XY = get_test_aucs(XY_results[, 1], XY_results[, -1]),
    CY_X = get_test_aucs(CY_X_results[, 1], CY_X_results[, -1]),
    lcd = get_lcd_aucs(labels_lcd, CX_results[, -1], XY_results[, -1], CY_X_results[, -1]),
    times = times_trans))
}

experiment_lcd_compare_tests_auc <- function(m = 200, dim_C = 2,
                     Ns = c(60, 80, 100,
                            seq(120, 180, by = 30),
                            seq(200, 500, by = 100),
                            600, 800, 1000, 1250, 1500),
                     simulation = NULL,
                     path = 'output/lcd_compare_tests/') {
  aucs <- experiment_lcd_compare_tests_roc(m = m, n = Ns[[1]], dim_C = dim_C, simulation = simulation, save_figures = FALSE)
  for (n in Ns[-1]) {
    res <- experiment_lcd_compare_tests_roc(m = m, n = n, dim_C = dim_C, simulation = simulation, save_figures = FALSE)
    for (type in names(aucs)) {
      for (test in names(aucs[[type]])) {
        aucs[[type]][[test]] <- c(aucs[[type]][[test]], res[[type]][[test]])
      }
    }
  }
  
  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
  
  save(aucs, file = sprintf("%saucs-%s.Rdata", path, timestamp))

  .ggsave(paste(path, 'CX,dim_C=', dim_C, sep = ""), .get_auc_plot(aucs$CX, Ns), 8, 8)
  .ggsave(paste(path, 'XY,dim_C=', dim_C, sep = ""), .get_auc_plot(aucs$XY, Ns), 8, 8)
  .ggsave(paste(path, 'CY_X,dim_C=', dim_C, sep = ""), .get_auc_plot(aucs$CY_X, Ns), 8, 8)
  .ggsave(paste(path, 'lcd,dim_C=', dim_C, sep = ""), .get_auc_plot(aucs$lcd, Ns), 8, 8)
  .ggsave(paste(path, 'time,dim_C=', dim_C, sep = ""), .plot_auc_times(aucs$times, Ns, lap = FALSE), 8, 8)
}
