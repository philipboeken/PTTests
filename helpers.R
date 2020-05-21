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

.hhist <- function(x) {
  ggplot(data.frame(x), aes(x)) + geom_histogram(bins=length(x)/10)
}

.gamplot <- function(x, y, binomial=FALSE) {
  if (binomial)
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x), method.args = list(family = "binomial"))
  else
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x))
}

.pplot3d <- function(X,Y,Z) {
  scene <- list(
    xaxis = list(title = substitute(X)),
    yaxis = list(title = substitute(Y)),
    zaxis = list(title = substitute(Z))
  )
  fig <- plot_ly(x=X, y=Y, z=Z, type="scatter3d", mode="markers", color=Z)
  fig <- fig %>% layout(scene=scene)
  fig
}

.pplotcs <- function(X, Y, Z, Z_min, Z_max) {
  XYZ <- cbind(X, Y, Z)
  data <- as.matrix(XYZ[Z_min < XYZ[,3] & XYZ[,3] < Z_max, c(1, 2)])
  pplot(data[,1], data[,2]) +
    labs(x = substitute(X), y = substitute(Y))
}

.scale_norm <- function(X) {
  return(pnorm(X, mean=mean(X), sd=sd(X)))
}

.scale_lin <- function(X) {
  if (ncol(as.matrix(X)) > 1) {
    X <- apply(X, 2, scale_lin)
  }
  return((X - min(X)) / (max(X) - min(X)))
}

pplot_roc <- function(labels, predictions, title=NULL, legend_pos=c(0.78, 0.25)) {
  predictions <- as.matrix(predictions)
  roc_data <- c()
  for (i in 1:ncol(predictions)) {
    pred <- prediction(predictions[,i], labels)
    res <- performance(pred, "tpr", "fpr")
    auc <- round(performance(pred, "auc")@y.values[[1]], 3)
    x <- res@x.values[[1]]
    y <- res@y.values[[1]]
    name <- colnames(predictions)[i]
    roc_data[[name]] <- list(data=data.frame(x=x, y=y), auc=auc, name=name)
  }
  
  plt <- ggplot() + 
    labs(x="False Positive Rate", y="True Positive Rate", title=title) +
    theme(legend.title = element_blank(),
          legend.position = legend_pos,
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- paste(roc$name, ' (', roc$auc, ')', sep="")
    plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}}))
  }
  
  return(plt)
}

pplot_roc_custom <- function(labels, ts_res, uci_res, ci_res, title=NULL, opt=1) {
  roc_data <- c()
  for (i in 1:ncol(ts_res)) {
    name <- colnames(ts_res)[i]
    bayes <- (name == 'bayes')
    res <- .lcd_performance(labels, ts_res[,i], uci_res[,i], ci_res[,i], bayes, opt)
    x <- res$fpr
    y <- res$tpr
    roc_data[[name]] <- list(data=data.frame(x=x, y=y), auc=res$auc, name=name)
  }
  
  plt <- ggplot() + 
    labs(x="False Positive Rate", y="True Positive Rate", title=title) +
    theme(legend.title = element_blank(),
          legend.position = c(0.78, 0.25),
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- paste(roc$name, ' (', roc$auc, ')', sep="")
    plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}}))
  }
  
  return(plt)
}

.lcd_performance <- function(labels, ts, uci, ci, bayes, opt=1) {
  
  alphas <- sort(c(Inf, ts, uci, ci), TRUE)
  
  false <- which(labels == 0)
  true <- which(labels == 1)
  
  n <- length(labels)
  a <- ifelse(bayes, 0.5, 1/(3*sqrt(n)))
  
  fp <- c()
  tp <- c()
  for (alpha in alphas) {
    if (opt) {
      idx <- intersect(intersect(which(ts >= alpha), which(uci >= alpha)), which(ci >= min(alpha, a)))
    } else {
      idx <- intersect(intersect(which(ts >= alpha), which(uci >= alpha)), which(ci >= a))
    }
    # idx <- intersect(intersect(which(ts >= alpha), which(uci >= alpha)), which(ci >= alpha))
    tp <- c(tp, length(intersect(idx, true)))
    fp <- c(fp, length(intersect(idx, false)))
  }
  
  tpr <- tp / length(true)
  fpr <- fp / length(false)
  
  auc <- .get_auc(tpr, fpr)
  # auc <- NaN
  
  return(list(tpr=c(tpr, 1), fpr=c(fpr, 1), auc=auc))
}

.get_auc <- function(tpr, fpr) {
  if (0 < min(tpr)) {
    tpr <- c(min(tpr)/2, min(tpr)/2, tpr)
    fpr <- c(0, min(fpr), fpr)
  }
  if (max(tpr) < 1) {
    tpr <- c(tpr, max(tpr) + (1-max(tpr))/2)
    fpr <- c(fpr, max(fpr) + (1-max(fpr))/2)
  }
  int <- round(sum(tpr[-1] * diff(fpr)), 3)
  return(int)
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




