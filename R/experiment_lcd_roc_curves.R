experiment_lcd_roc_curves <- function(m = 1000, n = 400, graph_probs = c(3 / 5, 1 / 5, 1 / 5),
                                       dim_C = 2, err_sd = 1 / 2, p_link = 4 / 5,
                                       interv_options = c(mean_shift, variance_shift, mixture),
                                       nonlin_options = c(linear, parabolic, sinusoidal),
                                       seed = 0,
                                       path = 'output/lcd_roc_curves/') {

  # Setup test
  ##############################################

  set.seed(seed)

  get_results <- function(dataset, test) {
    `%dopar%` <- foreach::`%dopar%`
    foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
      data <- dataset[[i]]

      CX <- test(data$C, data$X)
      XY <- test(data$X, data$Y)
      CY_X <- test(data$C, data$Y, data$X)

      data.frame(label_CX = data$label_CX,
                 label_XY = data$label_XY,
                 label_CY_X = data$label_CY_X,
                 label_lcd = data$label_lcd,
                 CX = CX, XY = XY, CY_X = CY_X
      )
    }
  }


  # Do test
  ##############################################

  doParallel::registerDoParallel()

  data <- lapply(1:m, function(i) get_data(graph_probs, n, dim_C, err_sd, p_link,
                                            interv_options, nonlin_options))

  results <- list(
    ppcor = get_results(data, .ppcor_wrapper),
    ppcor_b = get_results(data, .ppcor_b_wrapper),
    polyatree = get_results(data, .polyatree_wrapper)
  )


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

  CX_results <- .get_results_by_type(results, 'CX')
  XY_results <- .get_results_by_type(results, 'XY')
  CY_X_results <- .get_results_by_type(results, 'CY_X')

  labels_lcd <- factor(results$polyatree[, 'label_lcd'], ordered = TRUE, levels = c(1, 0))

  t0 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > \\min(\\alpha_0, 1-\\alpha))$')
  t1 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > \\alpha)$')
  t2 <- latex2exp::TeX('$(p_{CX} < \\alpha)$ and $(p_{XY} < \\alpha)$ and $(p_{CY|X} > 1-\\alpha)$')
  .plot0 <- .plot_roc_custom(labels_lcd, CX_results[, -1], XY_results[, -1], CY_X_results[, -1],
                             t0, option = 0, plot_point = FALSE)
  .plot1 <- .plot_roc_custom(labels_lcd, CX_results[, -1], XY_results[, -1], CY_X_results[, -1],
                             t1, option = 1, plot_point = FALSE)
  .plot2 <- .plot_roc_custom(labels_lcd, CX_results[, -1], XY_results[, -1], CY_X_results[, -1],
                             t2, option = 2, plot_point = FALSE)

  grid <- cowplot::plot_grid(.plot0, .plot1, .plot2, nrow = 1)

  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

  save(results, file = sprintf("%s%s.Rdata", path, timestamp))
  .ggsave(paste(path, timestamp, sep = ''), grid, 30, 10)
  .ggsave(paste(path, 'last', sep = ''), grid, 30, 10)
  .ggsave(paste(path, 'plot0', sep = ''), .plot0, 10, 10)
  .ggsave(paste(path, 'plot1', sep = ''), .plot1, 10, 10)
  .ggsave(paste(path, 'plot2', sep = ''), .plot2, 10, 10)
}
