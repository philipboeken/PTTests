.plot_colours <- list(polyatree = "#0072B2", polyatree_c = "#CC79A7",
                      ppcor = "#D55E00", spcor = "#56B4E9",
                      gcm = "#E69F00", rcot = "#009E73",
                      ccit = "#293352", ppcor_b = "#CC79A7",
                      mi_mixed = "#0072B2", lr_mixed = "#293352")
# https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/

.plot_roc <- function(labels, predictions, title = NULL, legend_pos = c(0.78, 0.275),
                      freq_default = 0.05, plot_point = TRUE) {
  predictions <- as.matrix(predictions)
  roc_data <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("Test", "fpr", "tpr")))
  dot_data <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("Test", "fpr", "tpr")))
  info <- list()
  for (i in 1:ncol(predictions)) {
    name <- colnames(predictions)[i]
    roc <- .get_roc(labels, - predictions[, i])
    bayes <- (name == 'polyatree' || name == 'polyatree_c' || name == 'ppcor_b')
    dot <- .get_roc_point(labels, predictions[, i], bayes, freq_default)
    info[[name]] <- ifelse(is.na(roc$auc), name, paste(name, ' (', roc$auc, ')', sep = ""))
    roc_data <- rbind(roc_data, cbind(rep(name, length(roc$x)), roc$x, roc$y))
    dot_data <- rbind(dot_data, cbind(name, dot$fpr, dot$tpr))
  }
  roc_data <- data.frame(Test = roc_data[, 1], fpr = as.numeric(roc_data[, 2]),
                         tpr = as.numeric(roc_data[, 3]))
  dot_data <- data.frame(Test = dot_data[, 1], fpr = as.numeric(dot_data[, 2]),
                         tpr = as.numeric(dot_data[, 3]))
  return(ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr, group = Test)) +
           ggplot2::geom_line(ggplot2::aes(colour = Test)) +
           ggplot2::scale_x_continuous("False Positive Rate", limits = c(0, 1)) +
           ggplot2::scale_y_continuous("True Positive Rate", limits = c(0, 1)) +
           ggplot2::theme(legend.title = ggplot2::element_text(size = 0),
                          legend.spacing.x = ggplot2::unit(0.2, 'cm'),
                          legend.position = legend_pos) +
           ggplot2::geom_point(data = dot_data, ggplot2::aes(colour = Test)) +
           ggplot2::scale_color_manual(values = unlist(.plot_colours), labels = unlist(info)))
}

.plot_roc_lcd <- function(labels, CX_results, XY_results, CY_X_results,
                          title = NULL, plot_point = TRUE, option = 0, legend_pos = c(0.78, 0.275)) {
  roc_data <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("Test", "fpr", "tpr")))
  dot_data <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("Test", "fpr", "tpr")))
  info <- list()
  for (i in 1:ncol(CX_results)) {
    name <- colnames(CX_results)[i]
    bayes <- (name == 'polyatree' || name == 'ppcor_b' || name == 'polyatree_c')
    roc <- .get_lcd_roc(labels, CX_results[, i], XY_results[, i], CY_X_results[, i], bayes, option)
    dot <- .get_lcd_roc_point(labels, CX_results[, i], XY_results[, i], CY_X_results[, i], bayes)
    info[[name]] <- name
    # info[[name]] <- ifelse(is.na(roc$auc), name, paste(name, ' (', roc$auc, ')', sep = ""))
    roc_data <- rbind(roc_data, cbind(rep(name, length(roc$fpr)), roc$fpr, roc$tpr))
    dot_data <- rbind(dot_data, cbind(name, dot$fpr, dot$tpr))
  }
  roc_data <- data.frame(Test = roc_data[, 1], fpr = as.numeric(roc_data[, 2]),
                         tpr = as.numeric(roc_data[, 3]))
  dot_data <- data.frame(Test = dot_data[, 1], fpr = as.numeric(dot_data[, 2]),
                         tpr = as.numeric(dot_data[, 3]))
  
  plt <- ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr, group = Test))
  if (option == 1) {
    plt <- plt + ggplot2::geom_path(ggplot2::aes(colour = Test))
  } else {
    plt <- plt + ggplot2::geom_line(ggplot2::aes(colour = Test))
  }
  
  if (plot_point) {
    plt <- plt + ggplot2::geom_point(data = dot_data, ggplot2::aes(colour = Test))
  }
  
  plt + ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous("False Positive Rate", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("True Positive Rate", limits = c(0, 1)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 0),
                   legend.spacing.x = ggplot2::unit(0.2, 'cm'),
                   legend.position = legend_pos,
                   plot.title = ggplot2::element_text(size = 11, hjust = 0.5)) +
    ggplot2::scale_color_manual(values = unlist(.plot_colours), labels = unlist(info))
}

.get_roc <- function(labels, predictions) {
  pred <- ROCR::prediction(predictions, labels)
  res <- ROCR::performance(pred, "tpr", "fpr")
  auc <- round(ROCR::performance(pred, "auc")@y.values[[1]], 3)
  list(x = res@x.values[[1]], y = res@y.values[[1]], auc = auc)
}

.get_roc_point <- function(labels, predictions, bayes, freq_default = 0.05) {
  true <- which(labels == 0)
  false <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, freq_default)
  idx <- which(predictions <= alpha)
  tpr <- length(intersect(idx, true)) / length(true)
  fpr <- length(intersect(idx, false)) / length(false)
  return(list(tpr = tpr, fpr = fpr))
}

.get_lcd_roc <- function(labels, CX_results, XY_results, CY_X_results, bayes, option = 0) {
  alphas <- sort(c(CX_results, XY_results, CY_X_results), TRUE)
  alphas <- alphas[alphas != 0]
  
  false <- which(labels == 0)
  true <- which(labels == 1)
  
  n <- length(labels)
  a0 <- ifelse(bayes, 0.5, 1 / (5 * sqrt(n)))
  
  fp <- c()
  tp <- c()
  for (alpha in rev(alphas)) {
    idx <- which(CX_results <= alpha)
    idx <- intersect(idx, which(XY_results <= alpha))
    if (option == 0) {
      idx <- intersect(idx, which(CY_X_results >= min(a0, 1 - alpha)))
    } else if (option == 1) {
      idx <- intersect(idx, which(CY_X_results >= alpha))
    } else if (option == 2) {
      idx <- intersect(idx, which(CY_X_results >= 1 - alpha))
    }
    
    tp <- c(tp, length(intersect(idx, true)))
    fp <- c(fp, length(intersect(idx, false)))
  }
  
  tpr <- tp / length(true)
  fpr <- fp / length(false)
  
  if (option == 1) {
    return(list(tpr = c(0, tpr), fpr = c(0, fpr), auc = NaN))
  } else {
    auc <- .get_auc(tpr, fpr)
    return(list(tpr = c(0, tpr, 1), fpr = c(0, fpr, 1), auc = auc))
  }
}

.get_auc_plot <- function(results, Ns) {
  auc_data <- lapply(results, function(result) unlist(result))
  auc_data$n <- Ns
  auc_data <- reshape::melt(data.frame(auc_data), id.vars = 'n', variable_name = "Test:")
  ggplot2::ggplot(auc_data, ggplot2::aes(n)) +
    ggplot2::scale_x_continuous(limits = c(min(Ns), max(Ns)), trans = scales::log10_trans()) +
    ggplot2::geom_line(ggplot2::aes(y = value, colour = `Test:`)) +
    ggplot2::scale_color_manual(values = unlist(.plot_colours)) +
    ggplot2::labs(x = "n", y = "AUC") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(legend.position = "none")
}

.get_lcd_roc_point <- function(labels, CX_results, XY_results, CY_X_results, bayes) {
  false <- which(labels == 0)
  true <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, 0.01)
  idx <- which(CX_results <= alpha)
  idx <- intersect(idx, which(XY_results <= alpha))
  idx <- intersect(idx, which(CY_X_results >= alpha))
  tpr <- length(intersect(idx, true)) / length(true)
  fpr <- length(intersect(idx, false)) / length(false)
  return(list(tpr = tpr, fpr = fpr))
}

.get_auc <- function(tpr, fpr) {
  if (0 < min(fpr)) {
    tpr <- c(min(tpr) / 2, min(tpr) / 2, tpr)
    fpr <- c(0, min(fpr), fpr)
  }
  if (max(fpr) < 1) {
    tpr <- c(tpr, 1 - (1 - max(tpr)) / 2, 1 - (1 - max(tpr)) / 2)
    fpr <- c(fpr, max(fpr), 1)
  }
  auc <- round(sum(tpr[-1] * diff(fpr)), 3)
  return(auc)
}

.ggsave <- function(name, grid, width, height) {
  ggplot2::ggsave(
    paste(name, ".pdf", sep = ""),
    plot = grid, scale = 1,
    width = width,
    height = height,
    units = "cm",
    limitsize = FALSE
  )
}

.graph_lines <- function(context, system, context_edges, system_edges, opts = "", opts_context = "", plot_context = FALSE) {
  opts <- ifelse(opts != "", sprintf(", %s", opts), "")
  opts_context <- ifelse(opts_context != "", sprintf(", %s", opts_context), "")
  output <- ""
  if (plot_context) {
    for (node in context) {
      output <- sprintf(
        "%s\n\"%s\"[label = \"%s\", shape = box%s];",
        output, node, node, opts_context)
    }
  }
  for (node in system) {
    output <- sprintf(
      "%s\n\"%s\"[label=\"%s\", shape=oval%s];",
      output, node, node, opts)
  }
  for (i in 1:nrow(system_edges)) {
    output <- sprintf(
      "%s\n\"%s\"->\"%s\"[arrowtail=\"none\", arrowhead=\"normal\"%s];",
      output, system_edges[i, 1], system_edges[i, 2], opts)
  }
  if (plot_context) {
    for (i in 1:nrow(context_edges)) {
      output <- sprintf(
        "%s\n\"%s\"->\"%s\"[arrowtail=\"none\", arrowhead=\"normal\", style=\"dashed\"%s];",
        output, context_edges[i, 1], context_edges[i, 2], opts_context)
    }
  }
  
  return(output)
}

.output_graph <- function(results, path, name, alpha1, alpha2, plot_context = FALSE) {
  options(warn = 2)
  strong <- dplyr::filter(results, CX <= alpha1$strong, XY <= alpha1$strong, CY_X >= alpha2$strong)
  strong <- strong[, c('C', 'X', 'Y')]
  substantial <- dplyr::filter(results, CX <= alpha1$substantial,
                               XY <= alpha1$substantial, CY_X >= alpha2$substantial)
  substantial <- substantial[, c('C', 'X', 'Y')]
  weak <- dplyr::filter(results, CX <= alpha1$weak, XY <= alpha1$weak, CY_X >= alpha2$weak)
  weak <- weak[, c('C', 'X', 'Y')]
  output <- "digraph G {margin=0;"
  
  context1 <- unique(as.character(strong[, 'C']))
  system1 <- unique(c(as.character(strong[, 'X']), as.character(strong[, 'Y'])))
  context_edges1 <- unique(strong[, c('C', 'X')])
  system_edges1 <- unique(strong[, c('X', 'Y')])
  output <- paste(output, .graph_lines(context1, system1, context_edges1, system_edges1,
                                       "color=\"#000000\"",
                                       "color=\"#c4c4c4\"",
                                       plot_context = plot_context), sep = "")
  
  diff <- dplyr::setdiff(substantial, strong)
  context2 <- dplyr::setdiff(unique(as.character(diff[, 'C'])), context1)
  system2 <- dplyr::setdiff(unique(c(as.character(diff[, 'X']), as.character(diff[, 'Y']))), system1)
  context_edges2 <- dplyr::setdiff(unique(diff[, c('C', 'X')]), context_edges1)
  system_edges2 <- dplyr::setdiff(unique(diff[, c('X', 'Y')]), system_edges1)
  output <- paste(output, .graph_lines(context2, system2, context_edges2, system_edges2,
                                       "color=\"#e30505\"",
                                       "color=\"#ffa1a1\"",
                                       plot_context = plot_context), sep = "")
  
  diff <- dplyr::setdiff(weak, substantial)
  context3 <- dplyr::setdiff(unique(as.character(diff[, 'C'])),
                             c(context1, context2))
  system3 <- dplyr::setdiff(unique(c(as.character(diff[, 'X']), as.character(diff[, 'Y']))),
                            c(system1, system2))
  context_edges3 <- dplyr::setdiff(unique(diff[, c('C', 'X')]),
                                   rbind(context_edges1, context_edges2))
  system_edges3 <- dplyr::setdiff(unique(diff[, c('X', 'Y')]),
                                  rbind(system_edges1, system_edges2))
  output <- paste(output, .graph_lines(context3, system3, context_edges3, system_edges3,
                                       "color=\"#0032fa\"",
                                       "color=\"#a8baff\"",
                                       plot_context = plot_context), sep = "")
  output <- paste(output, "}", sep = "\n")
  
  write(output, file = paste(path, name, ".dot", sep = ""))
}

.plot_roc_times <- function(times, title = NULL) {
  times$time <- sapply(times$time, function(x)(x < 1.1) * 1.1 + (x >= 1.1) * x)
  ggplot2::ggplot(data = times, ggplot2::aes(x = ensemble, y = time, fill = test)) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::labs(x = "Test ensemble", y = "Runtime (sec.)", title = title) +
    ggplot2::scale_fill_discrete(name = "",
                                 breaks = c("1_CX", "2_XY", "3_CY_X"),
                                 labels = c("Two-sample", "Indep.", "Cond. indep.")) +
    ggplot2::theme(legend.position = c(0.703, 0.915),
                   legend.direction = "horizontal") +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
                           labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))
}

.plot_auc_times <- function(results, Ns, save_legend = FALSE, lap = TRUE) {
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

.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

.scatterplot <- function(data) {
  plt <- ggplot2::ggplot() + ggplot2::geom_point(data = data, ggplot2::aes(x = Y, y = X, colour = C))
  if (length(unique(data$C)) == 2) {
    return(plt + ggplot2::scale_color_manual(labels = c("C = 0", "C = 1"), values = c("blue", "red")) +
             ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = c(0.87, 0.12)))
  } else {
    labels <- sapply(sort(unique(data$C)), function(c) paste("C = ", c, sep = ""))
    return(plt + ggplot2::scale_colour_manual(labels = labels, values = c("#00AFBB", "#E7B800", "#52854C")) +
             ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = c(0.50, 0.5)))
  }
}
