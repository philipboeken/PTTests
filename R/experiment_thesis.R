experiment_thesis <- function(m = 200, dim_C = 3, err_sd = 1 / 2,
                              Ns = c(25, 30, seq(40, 100, by = 15), seq(120, 180, by = 30), seq(200, 500, by = 100), 600, 800, 1000, 1250, 1500),
                              save_legend = TRUE) {
  discrete_v_continuous <- function(m, dim_C, err_sd, Ns,
                                    save_legend = FALSE, seed = 1, path = 'output/thesis/discrete/') {
    set.seed(seed)
    results_1 <- do_test(1, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 0, interv = c(mean_shift), link = c(linear),
                         var1 = 'X', var2 = 'C', var3 = NULL, save_legend = save_legend)
    results_2 <- do_test(2, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(mean_shift), link = c(linear),
                         var1 = 'X', var2 = 'C', var3 = NULL)
    results_3 <- do_test(3, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(variance_shift), link = c(linear),
                         var1 = 'X', var2 = 'C', var3 = NULL)
    results_4 <- do_test(4, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(mixture), link = c(linear),
                         var1 = 'X', var2 = 'C', var3 = NULL)
    
    save_results(results_1, results_2, results_3, results_4, path, 'discrete')
  }
  
  continuous_v_continuous <- function(m, dim_C, err_sd, Ns,
                                      seed = 0, path = 'output/thesis/continuous/') {
    set.seed(seed)
    results_1 <- do_test(1, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 0, interv = c(no_intervention), link = c(linear),
                         var1 = 'X', var2 = 'Y', var3 = NULL)
    results_2 <- do_test(2, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(no_intervention), link = c(linear),
                         var1 = 'X', var2 = 'Y', var3 = NULL)
    results_3 <- do_test(3, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(no_intervention), link = c(parabolic),
                         var1 = 'X', var2 = 'Y', var3 = NULL)
    results_4 <- do_test(4, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(no_intervention), link = c(sinusoidal),
                         var1 = 'X', var2 = 'Y', var3 = NULL)
    
    save_results(results_1, results_2, results_3, results_4, path, 'continuous')
  }
  
  discrete_v_continuous_ci <- function(m, dim_C, err_sd, Ns,
                                       save_legend = FALSE, seed = 0, path = 'output/thesis/discrete_ci/') {
    set.seed(seed)
    results_1 <- do_test(1, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(mean_shift), link = c(parabolic),
                         var1 = 'C', var2 = 'Y', var3 = 'X', save_legend = save_legend)
    results_2 <- do_test(2, m, dim_C, err_sd, Ns, graph_probs = c(0, 1, 0),
                         p_link = 1, interv = c(mean_shift), link = c(parabolic),
                         var1 = 'C', var2 = 'Y', var3 = 'X')
    results_3 <- do_test(3, m, dim_C, err_sd, Ns, graph_probs = c(1, 0, 0),
                         p_link = 1, interv = c(mean_shift), link = c(sinusoidal),
                         var1 = 'C', var2 = 'Y', var3 = 'X')
    results_4 <- do_test(4, m, dim_C, err_sd, Ns, graph_probs = c(0, 1, 0),
                         p_link = 1, interv = c(mean_shift), link = c(sinusoidal),
                         var1 = 'C', var2 = 'Y', var3 = 'X')
    
    save_results(results_1, results_2, results_3, results_4, path, 'discrete_ci')
  }
  
  save_results <- function(results_1, results_2, results_3, results_4, path, name) {
    .ggsave(paste(path, '1_data', sep = ''), results_1$data_plot, 7, 7)
    .ggsave(paste(path, '1_p_h0', sep = ''), results_1$pt_plot, 7, 7)
    .ggsave(paste(path, '1_pval', sep = ''), results_1$pval_plot, 7, 7)
    .ggsave(paste(path, '1_time', sep = ''), results_1$time_plot, 7, 7)
    
    .ggsave(paste(path, '2_data', sep = ''), results_2$data_plot, 7, 7)
    .ggsave(paste(path, '2_p_h0', sep = ''), results_2$pt_plot, 7, 7)
    .ggsave(paste(path, '2_pval', sep = ''), results_2$pval_plot, 7, 7)
    .ggsave(paste(path, '2_time', sep = ''), results_2$time_plot, 7, 7)
    
    .ggsave(paste(path, '3_data', sep = ''), results_3$data_plot, 7, 7)
    .ggsave(paste(path, '3_p_h0', sep = ''), results_3$pt_plot, 7, 7)
    .ggsave(paste(path, '3_pval', sep = ''), results_3$pval_plot, 7, 7)
    .ggsave(paste(path, '3_time', sep = ''), results_3$time_plot, 7, 7)
    
    .ggsave(paste(path, '4_data', sep = ''), results_4$data_plot, 7, 7)
    .ggsave(paste(path, '4_p_h0', sep = ''), results_4$pt_plot, 7, 7)
    .ggsave(paste(path, '4_pval', sep = ''), results_4$pval_plot, 7, 7)
    .ggsave(paste(path, '4_time', sep = ''), results_4$time_plot, 7, 7)
    
    grid <- cowplot::plot_grid(results_1$data_plot, results_2$data_plot,
                               results_3$data_plot, results_4$data_plot,
                               results_1$pt_plot, results_2$pt_plot,
                               results_3$pt_plot, results_4$pt_plot,
                               results_1$pval_plot, results_2$pval_plot,
                               results_3$pval_plot, results_4$pval_plot,
                               results_1$time_plot, results_2$time_plot,
                               results_3$time_plot, results_4$time_plot,
                               nrow = 4)
    
    .ggsave(paste(path, name, sep = ''), grid, 30, 30)
  }
  
  do_test <- function(testnr, m, dim_C, err_sd, Ns,
                      graph_probs, p_link, interv, link, var1, var2, var3 = NULL, save_legend = FALSE) {
    
    doParallel::registerDoParallel()
    
    results <- list(
      ppcor = list(test = .ppcor_wrapper, quantiles = c(), time = c()),
      spcor = list(test = .spcor_wrapper, quantiles = c(), time = c()),
      gcm = list(test = .gcm_wrapper, quantiles = c(), time = c()),
      rcot = list(test = .rcot_wrapper, quantiles = c(), time = c()),
      ccit = list(test = .ccit_wrapper, quantiles = c(), time = c()),
      polyatree = list(test = .pt_wrapper, quantiles = c(), time = c())
    )
    
    for (n in Ns) {
      cat(format(Sys.time(), "%X"), "Running test", testnr, "with n =", n, "\n")
      dataset <- get_data_transformed(m, graph_probs, n, dim_C, err_sd, p_link, interv, link, var1, var2, var3)
      for (test in names(results)) {
        `%dopar%` <- foreach::`%dopar%`
        result <- foreach::foreach(i = 1:length(dataset), .combine = rbind) %dopar% {
          data <- dataset[[i]]
          
          start_time <- Sys.time()
          result <- results[[test]]$test(data$X, data$Y, data$Z)
          end_time <- Sys.time()
          
          return(data.frame(values = result, time = end_time - start_time))
        }
        results[[test]]$quantiles <- rbind(results[[test]]$quantiles,
                                           quantile(result$values, c(0.05, 0.25, 0.5, 0.75, 0.95)))
        results[[test]]$time <- c(results[[test]]$time, sum(result$time))
      }
    }
    
    data <- get_data_transformed(1, graph_probs, 200, dim_C, err_sd,
                                 p_link, interv, link, var1, var2, var3)[[1]]
    if (is.null(var3)) {
      data_plot <- get_data_plot(data$X, data$Y, var1, var2)
    } else {
      data <- data.frame(C = as.factor(data$X), Y = data$Y, X = data$Z)
      data_plot <- .scatterplot(data)
      if (save_legend) {
        .ggsave('output/thesis/legend_data', cowplot::get_legend(
          data_plot + ggplot2::theme(legend.direction = "horizontal") +
            ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))), 10, 1.2)
      }
      data_plot <- data_plot + ggplot2::theme(legend.position = "none")
    }
    pt_plot <- get_pt_plot(results, Ns)
    time_plot <- get_time_plot(results, Ns, save_legend)
    pval_plot <- get_pval_plot(results, Ns)
    
    return(list(data_plot = data_plot, pt_plot = pt_plot,
                time_plot = time_plot, pval_plot = pval_plot))
  }
  
  get_data_transformed <- function(m, graph_probs, n, dim_C, err_sd, p_link, interv, link, var1, var2, var3) {
    lapply(1:m, function(i) {
      data <- get_data(graph_probs, n, dim_C, err_sd, p_link, interv, link)
      res <- list(X = data[[var1]], Y = data[[var2]])
      if (!is.null(var3)) {
        res$Z <- data[[var3]]
      }
      return(res)
    })
  }
  
  get_data_plot <- function(X, Y, var1, var2) {
    ggplot2::ggplot(data = data.frame(X = X, Y = Y)) +
      ggplot2::geom_point(ggplot2::aes(x = X, y = Y), color = "#00AFBB") +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = c(0.87, 0.12)) +
      ggplot2::labs(x = var1, y = var2)
  }
  
  get_pt_plot <- function(results, Ns) {
    ggplot2::ggplot(data.frame(cbind(results$polyatree$quantiles, Ns)), ggplot2::aes(Ns)) +
      ggplot2::geom_line(ggplot2::aes(y = X50.), colour = "#0072B2") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = X5., ymax = X95.), fill = "#0072B2", alpha = 0.3) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = X25., ymax = X75.), fill = "#0072B2", alpha = 0.5) +
      ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
      ggplot2::ylim(0, 1) +
      ggplot2::labs(x = "n", y = latex2exp::TeX("P(H_0 | data)")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = unlist(.plot_colours))
  }
  
  get_pval_plot <- function(results, Ns) {
    pval_data <- lapply(results, function(result) result$quantiles[, '50%'])
    pval_data$polyatree <- NULL
    pval_data$n <- Ns
    pval_data <- reshape::melt(data.frame(pval_data), id.vars = 'n')
    ggplot2::ggplot(pval_data, ggplot2::aes(n)) +
      ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
      ggplot2::geom_line(ggplot2::aes(y = value, colour = variable)) +
      ggplot2::ylim(0, 1) +
      ggplot2::labs(x = "n", y = latex2exp::TeX("$p$-value")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = unlist(.plot_colours))
  }
  
  get_time_plot <- function(results, Ns, save_legend = FALSE, lap = TRUE) {
    if (lap) {
      time_data <- lapply(results, function(result) result$time)
    } else {
      time_data <- results
    }
    time_data$n <- Ns
    time_data <- reshape::melt(data.frame(time_data), id.vars = 'n', variable_name = "Test:")
    plt <- ggplot2::ggplot(time_data, ggplot2::aes(n)) +
      ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
      ggplot2::scale_y_continuous(trans = scales::log10_trans()) +
      ggplot2::geom_line(ggplot2::aes(y = value, colour = `Test:`)) +
      ggplot2::scale_color_manual(values = unlist(.plot_colours)) +
      ggplot2::labs(x = "n", y = "time (s)")
    if (save_legend) {
      .ggsave('output/thesis/legend', cowplot::get_legend(
        plt + ggplot2::theme(legend.direction = "horizontal") +
          ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))), 12, 1.2)
    }
    plt + ggplot2::theme(legend.position = "none")
    
  }
  
  get_auc_plot <- function(results, Ns) {
    auc_data <- lapply(results, function(result) unlist(result))
    auc_data$n <- Ns
    auc_data <- reshape::melt(data.frame(auc_data), id.vars = 'n', variable_name = "Test:")
    ggplot2::ggplot(auc_data, ggplot2::aes(n)) +
      ggplot2::scale_x_continuous(limits = c(25, max(Ns)), trans = scales::log10_trans()) +
      ggplot2::geom_line(ggplot2::aes(y = value, colour = `Test:`)) +
      ggplot2::scale_color_manual(values = unlist(.plot_colours)) +
      ggplot2::labs(x = "n", y = "AUC") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme(legend.position = "none")
  }
  
  lcd_aucs <- function(m, dim_C, Ns, path = 'output/thesis/compare_tests_auc/') {
    aucs <- experiment_lcd_compare_tests(m = m, n = Ns[[1]], dim_C = dim_C, save_figures = FALSE)
    for (n in Ns[-1]) {
      res <- experiment_lcd_compare_tests(m = m, n = n, dim_C = dim_C, save_figures = FALSE)
      for (type in names(aucs)) {
        for (test in names(aucs[[type]])) {
          aucs[[type]][[test]] <- c(aucs[[type]][[test]], res[[type]][[test]])
        }
      }
    }
    
    .ggsave(paste(path, 'CX,dim_C=', dim_C, sep = ""), get_auc_plot(aucs$CX, Ns), 7, 7)
    .ggsave(paste(path, 'XY,dim_C=', dim_C, sep = ""), get_auc_plot(aucs$XY, Ns), 7, 7)
    .ggsave(paste(path, 'CY_X,dim_C=', dim_C, sep = ""), get_auc_plot(aucs$CY_X, Ns), 7, 7)
    .ggsave(paste(path, 'lcd,dim_C=', dim_C, sep = ""), get_auc_plot(aucs$lcd, Ns), 7, 7)
    .ggsave(paste(path, 'time,dim_C=', dim_C, sep = ""), get_time_plot(aucs$times, Ns, lap = FALSE), 7, 7)
  }
  
  library(RCIT)
  
  cat(format(Sys.time(), "%X"), "Discrete vs continuous\n")
  discrete_v_continuous(m, dim_C, err_sd, Ns, save_legend = save_legend)
  
  cat(format(Sys.time(), "%X"), "Continuous vs continuous\n")
  continuous_v_continuous(m, dim_C, err_sd, Ns)
  
  cat(format(Sys.time(), "%X"), "Discrete vs continuous CI\n")
  discrete_v_continuous_ci(m, dim_C, err_sd, Ns, save_legend = save_legend)
  
  cat(format(Sys.time(), "%X"), "LCD ROC curves\n")
  experiment_lcd_roc_curves(path = 'output/thesis/compare_roc_curves/')
  
  cat(format(Sys.time(), "%X"), "LCD compare tests\n")
  experiment_lcd_compare_tests(dim_C = dim_C, path = 'output/thesis/compare_tests_roc/', pt_continuous = TRUE)
  
  cat(format(Sys.time(), "%X"), "LCD AUC scores, d=2\n")
  lcd_aucs(m, dim_C = 2, Ns)
  
  cat(format(Sys.time(), "%X"), "LCD AUC scores, d=4\n")
  lcd_aucs(m, dim_C = 4, Ns)
  
  cat(format(Sys.time(), "%X"), "LCD AUC scores, d=8\n")
  lcd_aucs(m, dim_C = 8, Ns)
  
  # cat(format(Sys.time(), "%X"), "Sachs\n")
  # experiment_sachs_lcd()
}

