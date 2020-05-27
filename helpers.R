suppressMessages(library(ggplot2))
require(scales)
suppressMessages(library(plotly))
suppressMessages(library(ROCR))
suppressWarnings(library(dplyr))

.pplot <- function(x, y, count=FALSE) {
  if (count) {
    ggplot(data.frame(y,x), aes(x, y)) + geom_count()
  } else {
    result <- bayes.UCItest(x, y)
    cat("\nBF(H0, H1):\t", result$bf, "\nP(H0 | XY):\t", result$p_H0, "\n")
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() +
      labs(x = substitute(x), y = substitute(y))
  }
}

pplot_roc <- function(labels, predictions, title=NULL, legend_pos=c(0.78, 0.25), freq_default=0.01) {
  predictions <- as.matrix(predictions)
  roc_data <- c()
  for (i in 1:ncol(predictions)) {
    name <- colnames(predictions)[i]
    pred <- prediction(-predictions[,i], labels)
    res <- performance(pred, "tpr", "fpr")
    auc <- round(performance(pred, "auc")@y.values[[1]], 3)
    x <- res@x.values[[1]]
    y <- res@y.values[[1]]
    bayes <- (name == 'bayes')
    dot <- .roc_dot(labels, predictions[,i], bayes, freq_default)
    info <- paste(name, ' (', auc, ')', sep="")
    roc_data[[name]] <- list(data=data.frame(x=x, y=y), 
                             point=data.frame(x=dot$fpr, y=dot$tpr), 
                             info=info)
  }
  
  plt <- ggplot() + 
    labs(x="False Positive Rate", y="True Positive Rate", title=title) +
    theme(legend.title = element_blank(),
          legend.position = legend_pos,
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- roc$info
    plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}})) + 
      geom_point(data=roc$point, aes(x, y, colour={{c}}))
  }
  
  return(plt)
}

pplot_roc_custom <- function(labels, ts_res, uci_res, ci_res, 
                             title=NULL, plot_point=TRUE, option=0) {
  roc_data <- c()
  for (i in 1:ncol(ts_res)) {
    name <- colnames(ts_res)[i]
    bayes <- (name == 'polyatree' || name == 'ppcor_b')
    roc <- .lcd_roc(labels, ts_res[,i], uci_res[,i], ci_res[,i], bayes, option)
    dot <- .lcd_roc_dot(labels, ts_res[,i], uci_res[,i], ci_res[,i], bayes)
    info <- paste(name, ' (', roc$auc, ')', sep="")
    roc_data[[name]] <- list(data=data.frame(x=roc$fpr, y=roc$tpr), 
                             point=data.frame(x=dot$fpr, y=dot$tpr), 
                             info=info)
  }
  
  plt <- ggplot() +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(x="False Positive Rate", y="True Positive Rate", title=title) +
    theme(legend.title = element_blank(),
          legend.position = c(0.78, 0.25),
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- roc$info
    
    if (option ==1) {
      plt <- plt + geom_point(data=roc$data, aes(x, y, colour={{c}}), size=0.1)
    } else {
      plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}}))
    }
    
    if (plot_point) {
      plt <- plt + geom_point(data=roc$point, aes(x, y, colour={{c}}))
    }
  }
  
  return(plt)
}

.lcd_roc <- function(labels, ts, uci, ci, bayes, option=0) {
  alphas <- sort(c(ts, uci, ci), TRUE)
  alphas <- alphas[alphas != 0]
  
  false <- which(labels == 0)
  true <- which(labels == 1)
  
  n <- length(labels)
  a0 <- ifelse(bayes, 0.5, 1/(5*sqrt(n)))
  
  fp <- c()
  tp <- c()
  for (alpha in rev(alphas)) {
    idx <- which(ts <= alpha)
    idx <- intersect(idx, which(uci <= alpha))
    if (option == 0) {
      idx <- intersect(idx, which(ci >= min(a0, 1-alpha)))
    } else if (option == 1) {
      idx <- intersect(idx, which(ci >= alpha))
    } else if (option == 2) {
      idx <- intersect(idx, which(ci >= 1-alpha))
    } else if (option == 3) {
      idx <- intersect(idx, which(ci >= a0))
    } else if (option == 4) {
      idx <- which(ts <= a0)
      idx <- intersect(idx, which(uci <= a0))
      idx <- intersect(idx, which(ci >= 1-alpha))
    }
    
    tp <- c(tp, length(intersect(idx, true)))
    fp <- c(fp, length(intersect(idx, false)))
  }
  
  tpr <- tp / length(true)
  fpr <- fp / length(false)
  
  if (option == 1) {
    return(list(tpr=c(0, tpr), fpr=c(0, fpr), auc="-"))
  } else {
    auc <- .get_auc(tpr, fpr)
    return(list(tpr=c(0, tpr, 1), fpr=c(0, fpr, 1), auc=auc))
  }
}

.roc_dot <- function(labels, predictions, bayes, freq_default=0.01) {
  true <- which(labels == 0)
  false <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, freq_default)
  idx <- which(predictions <= alpha)
  tpr <- length(intersect(idx, true)) / length(true)
  fpr <- length(intersect(idx, false)) / length(false)
  return(list(tpr=tpr, fpr=fpr))
}

.lcd_roc_dot <- function(labels, ts, uci, ci, bayes) {
  false <- which(labels == 0)
  true <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, 0.01)
  idx <- which(ts <= alpha)
  idx <- intersect(idx, which(uci <= alpha))
  idx <- intersect(idx, which(ci >= alpha))
  tpr <- length(intersect(idx, true)) / length(true)
  fpr <- length(intersect(idx, false)) / length(false)
  return(list(tpr=tpr, fpr=fpr))
}

.get_auc <- function(tpr, fpr) {
  if (0 < min(fpr)) {
    tpr <- c(min(tpr)/2, min(tpr)/2, tpr)
    fpr <- c(0, min(fpr), fpr)
  }
  if (max(fpr) < 1) {
    tpr <- c(tpr, 1 - (1-max(tpr))/2, 1 - (1-max(tpr))/2)
    fpr <- c(fpr, max(fpr), 1)
  }
  auc <- round(sum(tpr[-1] * diff(fpr)), 3)
  return(auc)
}

.ggsave <- function (name, grid, width, height) {
  ggsave(
    paste(name, ".pdf", sep=""),
    plot = grid, scale = 1,
    width = width,
    height = height,
    units = "cm", dpi = 300, limitsize = FALSE
  )
}

.graph_lines <- function(context, system, context_edges, system_edges, opts="", opts_context="") {
  opts <- ifelse(opts != "", sprintf(", %s", opts), "")
  opts_context <- ifelse(opts_context != "", sprintf(", %s", opts_context), "")
  output <- ""
  for (node in context) {
    output <- sprintf(
      "%s\n\"%s\"[label=\"%s\", shape=box%s];", 
      output, node, node, opts_context)
  }
  for (node in system) {
    output <- sprintf(
      "%s\n\"%s\"[label=\"%s\", shape=oval%s];", 
      output, node, node, opts)
  }
  for (i in 1:nrow(system_edges)) {
    output <- sprintf(
      "%s\n\"%s\"->\"%s\"[arrowtail=\"none\", arrowhead=\"normal\"%s];",
      output,
      system_edges[i, 1],
      system_edges[i, 2],
      opts)
  }
  for (i in 1:nrow(context_edges)) {
    output <- sprintf(
      "%s\n\"%s\"->\"%s\"[arrowtail=\"none\", arrowhead=\"normal\"%s];",
      output,
      context_edges[i, 1],
      context_edges[i, 2],
      opts_context)
  }
  
  return(output)
}

.output_graph <- function(lcd_triples, path, name) {
  context <- unique(as.character(lcd_triples[,'C']))
  system <- unique(c(as.character(lcd_triples[,'X']), as.character(lcd_triples[,'Y'])))
  context_edges <- unique(lcd_triples[,c('C', 'X')])
  system_edges <- unique(lcd_triples[,c('X', 'Y')])
  output <- "digraph G {"
  output <- paste(output, .graph_lines(context, system, context_edges, system_edges), sep="")
  output <- paste(output, "}", sep="\n")
  
  write(output, file=paste(path, name, ".dot", sep=""))
}

.output_graph_levels <- function(strong, substantial, weak, path, name) {
  strong <- strong[,c('C', 'X', 'Y')]
  substantial <- substantial[,c('C', 'X', 'Y')]
  weak <- weak[,c('C', 'X', 'Y')]
  output <- "digraph G {"
  
  context1 <- unique(as.character(strong[,'C']))
  system1 <- unique(c(as.character(strong[,'X']), as.character(strong[,'Y'])))
  context_edges1 <- unique(strong[,c('C', 'X')])
  system_edges1 <- unique(strong[,c('X', 'Y')])
  output <- paste(output, .graph_lines(context1, system1, context_edges1, system_edges1, 
                                       "color=\"#000000\"",
                                       "color=\"#c4c4c4\""), sep="")
  
  diff <- dplyr::setdiff(substantial, strong)
  context2 <- dplyr::setdiff(unique(as.character(diff[,'C'])), context1)
  system2 <- dplyr::setdiff(unique(c(as.character(diff[,'X']), as.character(diff[,'Y']))), system1)
  context_edges2 <- dplyr::setdiff(unique(diff[,c('C', 'X')]), context_edges1)
  system_edges2 <- dplyr::setdiff(unique(diff[,c('X', 'Y')]), system_edges1)
  output <- paste(output, .graph_lines(context2, system2, context_edges2, system_edges2, 
                                       "color=\"#e30505\"",
                                       "color=\"#ffa1a1\""), sep="")
  
  diff <- dplyr::setdiff(weak, substantial)
  context3 <- dplyr::setdiff(unique(as.character(diff[,'C'])), 
                             c(context1, context2))
  system3 <- dplyr::setdiff(unique(c(as.character(diff[,'X']), as.character(diff[,'Y']))), 
                            c(system1, system2))
  context_edges3 <- dplyr::setdiff(unique(diff[,c('C', 'X')]), 
                                   rbind(context_edges1, context_edges2))
  system_edges3 <- dplyr::setdiff(unique(diff[,c('X', 'Y')]), 
                                  rbind(system_edges1, system_edges2))
  output <- paste(output, .graph_lines(context3, system3, context_edges3, system_edges3, 
                                       "color=\"#0032fa\"",
                                       "color=\"#a8baff\""), sep="")
  output <- paste(output, "}", sep="\n")
  
  write(output, file=paste(path, name, ".dot", sep=""))
}

.output_sachs_plots <- function(all_data_pooled, results, filename) {  
  .plot_lcd_triple <- function(C, X, Y, title) {
    idx <- which(all_data_pooled[,C] == 1 | all_data_pooled[,'experiment'] == 1)
    data <- as.data.frame(all_data_pooled[idx, c(C, X, Y)])
    colnames(data) <- c('C', 'X', 'Y')
    data[,'C'] <- as.factor(data[,'C'])
    return(ggplot() + geom_point(data=data, aes(x=X, y=Y, colour=C)) + 
             labs(x="X", y="Y", title=title) +
             theme(legend.title = element_blank(),
                   legend.position = c(0.87, 0.12)) +
             scale_color_manual(labels = c("C=0", "C=1"), values = c("blue", "red")))
  }
  
  plots <- list()
  for (i in 1:nrow(results)) {
    C <- as.character(results[i,'C'])
    X <- as.character(results[i,'X'])
    Y <- as.character(results[i,'Y'])
    plots[[i]] <- .plot_lcd_triple(
      C, X, Y, sprintf("%s -> %s -> %s\nCX: %.3f, XY:  %.3f, CY|X:  %.3f",
                       C, X, Y, results[i,'CX'], results[i,'XY'], results[i,'CY_X']))
  }
  
  grid <- plot_grid(plotlist=as.list(plots), ncol=5, nrow=ceiling(nrow(results)/5))
  .ggsave(filename, grid, 50, ceiling(nrow(results)/5)*12)
}

plot_times <- function(times, title=NULL) {
  return(ggplot(data=times, aes(x=ensemble, y=time, fill=test)) +
           geom_bar(stat="identity", position=position_dodge()) +
           labs(x="Test ensemble", y="Runtime (sec.)", title=title) +
           scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))
  )
}
