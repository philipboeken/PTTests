suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(ROCR))

.pplot <- function(x, y, count=FALSE) {
  if (count) {
    ggplot(data.frame(y,x), aes(x, y)) + geom_count()
  } else {
    result <- bayes.UCItest(x, y)
    cat("\nBF(H0, H1):\t", result$bf, "\nP(H1 | XY):\t", result$p_H1, "\n")
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() +
      labs(x = substitute(x), y = substitute(y))
  }
}

pplot_roc <- function(labels, predictions, title=NULL, legend_pos=c(0.78, 0.25), freq_default=0.01) {
  predictions <- as.matrix(predictions)
  roc_data <- c()
  for (i in 1:ncol(predictions)) {
    name <- colnames(predictions)[i]
    pred <- prediction(predictions[,i], labels)
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

pplot_roc_custom <- function(labels, ts_res, uci_res, ci_res, title=NULL, opt=1) {
  roc_data <- c()
  for (i in 1:ncol(ts_res)) {
    name <- colnames(ts_res)[i]
    bayes <- (name == 'bayes')
    roc <- .lcd_roc(labels, ts_res[,i], uci_res[,i], ci_res[,i], bayes, opt)
    dot <- .lcd_roc_dot(labels, ts_res[,i], uci_res[,i], ci_res[,i], bayes)
    info <- paste(name, ' (', roc$auc, ')', sep="")
    roc_data[[name]] <- list(data=data.frame(x=roc$fpr, y=roc$tpr), 
                             point=data.frame(x=dot$fpr, y=dot$tpr), 
                             info=info)
  }
  
  plt <- ggplot() + 
    labs(x="False Positive Rate", y="True Positive Rate", title=title) +
    theme(legend.title = element_blank(),
          legend.position = c(0.78, 0.25),
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- roc$info
    plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}})) + 
      geom_point(data=roc$point, aes(x, y, colour={{c}}))
  }
  
  return(plt)
}

.lcd_roc <- function(labels, ts, uci, ci, bayes, opt=1) {
  alphas <- sort(c(ts, uci, ci), TRUE)
  alphas <- alphas[alphas != 0]
  
  false <- which(labels == 0)
  true <- which(labels == 1)
  
  n <- length(labels)
  a <- ifelse(bayes, 0.5, 1/(5*sqrt(n)))
  
  fp <- c()
  tp <- c()
  for (alpha in rev(alphas)) {
    
    if (opt == 1) {
      idx <- intersect(intersect(which(1-ts <= alpha), which(1-uci <= alpha)), which(ci >= min(a, 1-alpha)))
    } else if (opt == 2) {
      idx <- intersect(intersect(which(ts >= alpha), which(uci >= alpha)), which(ci >= a))
    } else {
      idx <- intersect(intersect(which(ts >= alpha), which(uci >= alpha)), which(ci >= alpha))
    }
    
    tp <- c(tp, length(intersect(idx, true)))
    fp <- c(fp, length(intersect(idx, false)))
  }
  
  tpr <- tp / length(true)
  fpr <- fp / length(false)
  
  auc <- .get_auc(tpr, fpr)
  
  return(list(tpr=c(0, tpr, 1), fpr=c(0, fpr, 1), auc=auc))
}

.roc_dot <- function(labels, predictions, bayes, freq_default=0.01) {
  false <- which(labels == 0)
  true <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, 1-freq_default)
  idx <- which(predictions >= alpha)
  tpr <- length(intersect(idx, true)) / length(true)
  fpr <- length(intersect(idx, false)) / length(false)
  return(list(tpr=tpr, fpr=fpr))
}

.lcd_roc_dot <- function(labels, ts, uci, ci, bayes) {
  false <- which(labels == 0)
  true <- which(labels == 1)
  n <- length(labels)
  alpha <- ifelse(bayes, 0.5, 0.01)
  idx <- intersect(intersect(which(1-ts <= alpha), which(1-uci <= alpha)), which(ci >= alpha))
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
    units = "cm", dpi = 300, limitsize = TRUE
  )
}




