experiment_thesis <- function(m = 60, dim_C = 4, err_sd = 1 / 2,
                              Ns = c(seq(40, 100, by = 20), seq(200, 500, by = 100), 600, 800, 1000),
                              save_legend = FALSE) {
  discrete_v_continuous <- function(m = 100, dim_C = 4, err_sd = 1 / 2,
                                    Ns = c(10, seq(20, 100, by = 20), seq(200, 500, by = 100), 600, 800, 1000),
                                    save_legend = FALSE, seed = 0, path = 'output/thesis/discrete/') {
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

    .ggsave(paste(path, '1_data', sep = ''), results_1$data_plot, 7, 7)
    .ggsave(paste(path, '1_p_h0', sep = ''), results_1$polyatree_plot, 7, 7)
    .ggsave(paste(path, '1_pval', sep = ''), results_1$pval_plot, 7, 7)
    # .ggsave(paste(path, '1_time', sep = ''), results_1$time_plot, 7, 7)

    .ggsave(paste(path, '2_data', sep = ''), results_2$data_plot, 7, 7)
    .ggsave(paste(path, '2_p_h0', sep = ''), results_2$polyatree_plot, 7, 7)
    .ggsave(paste(path, '2_pval', sep = ''), results_2$pval_plot, 7, 7)
    # .ggsave(paste(path, '2_time', sep = ''), results_2$time_plot, 7, 7)

    .ggsave(paste(path, '3_data', sep = ''), results_3$data_plot, 7, 7)
    .ggsave(paste(path, '3_p_h0', sep = ''), results_3$polyatree_plot, 7, 7)
    .ggsave(paste(path, '3_pval', sep = ''), results_3$pval_plot, 7, 7)
    # .ggsave(paste(path, '3_time', sep = ''), results_3$time_plot, 7, 7)

    .ggsave(paste(path, '4_data', sep = ''), results_4$data_plot, 7, 7)
    .ggsave(paste(path, '4_p_h0', sep = ''), results_4$polyatree_plot, 7, 7)
    .ggsave(paste(path, '4_pval', sep = ''), results_4$pval_plot, 7, 7)
    # .ggsave(paste(path, '4_time', sep = ''), results_4$time_plot, 7, 7)

    grid <- cowplot::plot_grid(results_1$data_plot, results_2$data_plot, results_3$data_plot, results_4$data_plot,
                               results_1$polyatree_plot, results_2$polyatree_plot, results_3$polyatree_plot, results_4$polyatree_plot,
                               results_1$pval_plot, results_2$pval_plot, results_3$pval_plot, results_4$pval_plot,
                               results_1$time_plot, results_2$time_plot, results_3$time_plot, results_4$time_plot,
                               nrow = 4)

    .ggsave(paste(path, 'discrete', sep = ''), grid, 30, 30)

  }

  continuous_v_continuous <- function(m = 100, dim_C = 4, err_sd = 1 / 2,
                                      Ns = c(10, seq(20, 100, by = 20), seq(200, 500, by = 100), 600, 800, 1000),
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

    .ggsave(paste(path, '1_data', sep = ''), results_1$data_plot, 7, 7)
    .ggsave(paste(path, '1_p_h0', sep = ''), results_1$polyatree_plot, 7, 7)
    .ggsave(paste(path, '1_pval', sep = ''), results_1$pval_plot, 7, 7)
    # .ggsave(paste(path, '1_time', sep = ''), results_1$time_plot, 7, 7)

    .ggsave(paste(path, '2_data', sep = ''), results_2$data_plot, 7, 7)
    .ggsave(paste(path, '2_p_h0', sep = ''), results_2$polyatree_plot, 7, 7)
    .ggsave(paste(path, '2_pval', sep = ''), results_2$pval_plot, 7, 7)
    # .ggsave(paste(path, '2_time', sep = ''), results_2$time_plot, 7, 7)

    .ggsave(paste(path, '3_data', sep = ''), results_3$data_plot, 7, 7)
    .ggsave(paste(path, '3_p_h0', sep = ''), results_3$polyatree_plot, 7, 7)
    .ggsave(paste(path, '3_pval', sep = ''), results_3$pval_plot, 7, 7)
    # .ggsave(paste(path, '3_time', sep = ''), results_3$time_plot, 7, 7)

    .ggsave(paste(path, '4_data', sep = ''), results_4$data_plot, 7, 7)
    .ggsave(paste(path, '4_p_h0', sep = ''), results_4$polyatree_plot, 7, 7)
    .ggsave(paste(path, '4_pval', sep = ''), results_4$pval_plot, 7, 7)
    # .ggsave(paste(path, '4_time', sep = ''), results_4$time_plot, 7, 7)

    grid <- cowplot::plot_grid(results_1$data_plot, results_2$data_plot, results_3$data_plot, results_4$data_plot,
                               results_1$polyatree_plot, results_2$polyatree_plot, results_3$polyatree_plot, results_4$polyatree_plot,
                               results_1$pval_plot, results_2$pval_plot, results_3$pval_plot, results_4$pval_plot,
                               results_1$time_plot, results_2$time_plot, results_3$time_plot, results_4$time_plot,
                               nrow = 4)

    timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

    .ggsave(paste(path, timestamp, sep = ''), grid, 40, 40)
    .ggsave(paste(path, 'continuous', sep = ''), grid, 40, 40)

  }

  discrete_v_continuous_ci <- function(m = 100, dim_C = 4, err_sd = 1 / 2,
                                       Ns = c(10, seq(20, 100, by = 20), seq(200, 500, by = 100), 600, 800, 1000),
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

    .ggsave(paste(path, '1_data', sep = ''), results_1$data_plot, 7, 7)
    .ggsave(paste(path, '1_p_h0', sep = ''), results_1$polyatree_plot, 7, 7)
    .ggsave(paste(path, '1_pval', sep = ''), results_1$pval_plot, 7, 7)
    # .ggsave(paste(path, '1_time', sep = ''), results_1$time_plot, 7, 7)

    .ggsave(paste(path, '2_data', sep = ''), results_2$data_plot, 7, 7)
    .ggsave(paste(path, '2_p_h0', sep = ''), results_2$polyatree_plot, 7, 7)
    .ggsave(paste(path, '2_pval', sep = ''), results_2$pval_plot, 7, 7)
    # .ggsave(paste(path, '2_time', sep = ''), results_2$time_plot, 7, 7)

    .ggsave(paste(path, '3_data', sep = ''), results_3$data_plot, 7, 7)
    .ggsave(paste(path, '3_p_h0', sep = ''), results_3$polyatree_plot, 7, 7)
    .ggsave(paste(path, '3_pval', sep = ''), results_3$pval_plot, 7, 7)
    # .ggsave(paste(path, '3_time', sep = ''), results_3$time_plot, 7, 7)

    .ggsave(paste(path, '4_data', sep = ''), results_4$data_plot, 7, 7)
    .ggsave(paste(path, '4_p_h0', sep = ''), results_4$polyatree_plot, 7, 7)
    .ggsave(paste(path, '4_pval', sep = ''), results_4$pval_plot, 7, 7)
    # .ggsave(paste(path, '4_time', sep = ''), results_4$time_plot, 7, 7)

    grid <- cowplot::plot_grid(results_1$data_plot, results_2$data_plot, results_3$data_plot, results_4$data_plot,
                               results_1$polyatree_plot, results_2$polyatree_plot, results_3$polyatree_plot, results_4$polyatree_plot,
                               results_1$pval_plot, results_2$pval_plot, results_3$pval_plot, results_4$pval_plot,
                               results_1$time_plot, results_2$time_plot, results_3$time_plot, results_4$time_plot,
                               nrow = 4)

    .ggsave(paste(path, 'discrete_ci', sep = ''), grid, 30, 30)

  }

  do_test <- function(testnr, m = 200, dim_C = 4, err_sd = 1 / 2, Ns = c(50, 100),
                      graph_probs, p_link, interv, link, var1, var2, var3 = NULL, save_legend = FALSE) {

    cores <- parallel::detectCores()
    cl <-parallel:: makeForkCluster(cores[1]-1)
    doParallel::registerDoParallel(cl)
    # doParallel::registerDoParallel()

    results <- list(
      ppcor = list(test = .ppcor_wrapper, quantiles = c(), time = c()),
      spcor = list(test = .spcor_wrapper, quantiles = c(), time = c()),
      gcm = list(test = .gcm_wrapper, quantiles = c(), time = c()),
      rcot = list(test = .rcot_wrapper, quantiles = c(), time = c()),
      ccit = list(test = .ccit_wrapper, quantiles = c(), time = c()),
      polyatree = list(test = .polyatree_wrapper, quantiles = c(), time = c())
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

    group_colors <- c(ppcor = "#E87D7E", spcor = "#B3A033", gcm = "#53B74C",
                      rcot = "#55BCC2", ccit = "#E46DDD", polyatree = "#6F9AF8")
    data <- get_data_transformed(1, graph_probs, 200, dim_C, err_sd, p_link, interv, link, var1, var2, var3)[[1]]
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
    polyatree_plot <- get_pt_plot(results, Ns, group_colors)
    time_plot <- get_time_plot(results, Ns, group_colors, save_legend)
    pval_plot <- get_pval_plot(results, Ns, group_colors)

    return(list(data_plot = data_plot, polyatree_plot = polyatree_plot, time_plot = time_plot, pval_plot = pval_plot))
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
    ggplot2::ggplot() +
      ggplot2::geom_point(data = data.frame(X = X, Y = Y), ggplot2::aes(x = X, y = Y)) +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = c(0.87, 0.12)) +
      ggplot2::labs(x = var1, y = var2)
  }

  get_pt_plot <- function(results, Ns, group_colors) {
    return(ggplot2::ggplot(data.frame(cbind(results$polyatree$quantiles, Ns)), ggplot2::aes(Ns)) +
             ggplot2::geom_line(ggplot2::aes(y = X50.), colour = "#6F9AF8") +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = X5., ymax = X95.), fill = "#6F9AF8", alpha = 0.3) +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = X25., ymax = X75.), fill = "#6F9AF8", alpha = 0.5) +
             ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
             ggplot2::ylim(0, 1) +
             ggplot2::labs(x = "n", y = latex2exp::TeX("P(H_0 | data)")) +
             ggplot2::theme(legend.position = "none") +
             ggplot2::scale_color_manual(values = group_colors)
    )
  }

  get_pval_plot <- function(results, Ns, group_colors) {
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
      ggplot2::scale_color_manual(values = group_colors)
  }

  get_time_plot <- function(results, Ns, group_colors, save_legend = FALSE) {
    time_data <- lapply(results, function(result) result$time)
    time_data$n <- Ns
    time_data <- reshape::melt(data.frame(time_data), id.vars = 'n', variable_name = "Test:")
    plt <- ggplot2::ggplot(time_data, ggplot2::aes(n)) +
      ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
      ggplot2::geom_line(ggplot2::aes(y = value, colour = `Test:`)) +
      ggplot2::scale_color_manual(values = group_colors) +
      ggplot2::labs(x = "n", y = "time (s)")
    if (save_legend) {
      .ggsave('output/thesis/legend', cowplot::get_legend(
        plt + ggplot2::theme(legend.direction = "horizontal") +
          ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))), 12, 1.2)
    }
    plt + ggplot2::theme(legend.position = "none")

  }

  discrete_v_continuous(m, dim_C, err_sd, Ns, save_legend=save_legend)
  continuous_v_continuous(m, dim_C, err_sd, Ns)
  discrete_v_continuous_ci(m, dim_C, err_sd, Ns, save_legend=save_legend)
}

# experiment_thesis(m=20, Ns = c(50, 60))

