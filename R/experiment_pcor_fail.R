experiment_pcor_fail <- function(m = 2000, n = 400, err_sd = 0.5,
                                 p_link = 0.8, p_C = 0.5,
                                 seed = 0,
                                 path = 'output/pcor_fail/') {

  # Setup test
  ##############################################

  set.seed(seed)

  get_data <- function(n, p_C, err_sd, p_link, link = 0) {
    # C -> X <- Y

    C <- rbinom(n, 1, p_C)
    Y <- rnorm(n)

    link_nonlin <- if (link < 0) 0 else if (link == 0) rbinom(1, 1, p_link) else 1
    X <- link_nonlin * Y ^ 2
    X <- X + err_sd * rnorm(n, 0, ifelse(sd(X) > 0, sd(X), 1 / err_sd))
    X <- C * X + (1 - C) * (X + 3)

    cond_indep <- as.numeric(!link_nonlin)

    return(list(C = C, Y = Y, X = X, label = cond_indep))
  }

  get_results <- function(m, n, p_C, err_sd, p_link, get_data) {
    `%dopar%` <- foreach::`%dopar%`
    result <- foreach::foreach(i = 1:m, .combine = rbind) %dopar% {
      data <- get_data(n, p_C, err_sd, p_link)
      return(data.frame(
        label = data$label,
        ppcor = .ppcor_wrapper(data$C, data$Y, data$X),
        polyatree = .polyatree_wrapper(data$C, data$Y, data$X)
      ))
    }
  }


  # Do test
  ##############################################

  doParallel::registerDoParallel()

  results <- get_results(m, n, p_C, err_sd, p_link, get_data)

  data_no_link <- get_data(n, p_C, err_sd, p_link, -1)
  data_linked <- get_data(n, p_C, err_sd, p_link, 1)


  # Process results
  ##############################################
  no_link <- data.frame(C = as.factor(data_no_link$C), Y = data_no_link$Y, X = data_no_link$X)
  linked <- data.frame(C = as.factor(data_linked$C), Y = data_linked$Y, X = data_linked$X)

  scat_plot_no_link <- .scatterplot(no_link)
  scat_plot_linked <- .scatterplot(linked)

  labels <- factor(results[, 1], ordered = TRUE, levels = c(1, 0))
  roc_plot <- .plot_roc(labels, results[, -1], NULL, c(0.8, 0.12), plot_point = FALSE)

  grid <- cowplot::plot_grid(scat_plot_no_link, scat_plot_linked, roc_plot, nrow = 1)

  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

  save(results, file = sprintf("%s%s.Rdata", path, timestamp))
  .ggsave(paste(path, 'ppcor_fail_no_link', sep = ''), scat_plot_no_link, 10, 10)
  .ggsave(paste(path, 'ppcor_fail_linked', sep = ''), scat_plot_linked, 10, 10)
  .ggsave(paste(path, 'ppcor_fail_roc', sep = ''), roc_plot, 10, 10)
  .ggsave(paste(path, timestamp, sep = ''), grid, 30, 10)
  .ggsave(paste(path, 'last', sep = ''), grid, 30, 10)
}
