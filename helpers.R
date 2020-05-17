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

pplot_roc <- function(labels, predictions, title=NULL) {
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
          legend.position = c(0.78, 0.25),
          plot.title = element_text(size=12, hjust=0.5))
  for (roc in roc_data) {
    c <- paste(roc$name, ' (', roc$auc, ')', sep="")
    plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}}))
  }
  
  return(plt)
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

# pplot_roc_custom <- function(labels, tests, title=NULL) {
#   tests <- as.matrix(tests)
#   pred <- .custom_lcd_prediction(labels, tests)
#   # res <- performance(pred, "tpr", "fpr")
#   # auc <- round(performance(pred, "auc")@y.values[[1]], 3)
#   x <- pred$fpr
#   y <- pred$tpr

#   plt <- ggplot() + 
#     labs(x="False Positive Rate", y="True Positive Rate", title=title) +
#     theme(legend.title = element_blank(),
#           legend.position = c(0.78, 0.25),
#           plot.title = element_text(size=12, hjust=0.5)) +
#     geom_line(data=data.frame(x=x,y=y), aes(x, y))

#   return(plt)
# }

# .custom_lcd_prediction <- function(labels, tests, method='freq') {
#   ts <- tests[,1]
#   uci <- tests[,2]
#   ci <- tests[,3]
#   score <- tests[,4]
#   alphas <- seq(0, 1, 1/(2*length(ts)-1))
#   fpr <- c(0, 1)
#   tpr <- c(0, 1)
#   false <- which(labels == 0)
#   true <- which(labels == 1)
#   for (alpha in alphas) {
#     # idx <- intersect(c(which(ts <= alpha), which(uci <= alpha)), which(ci >= 1-alpha))
#     cond <- which(score > alpha)
#     tpr <- c(tpr, length(intersect(cond, true)) / length(true))
#     fpr <- c(fpr, length(intersect(cond, false)) / length(ts))
#   }
#   return(list(tpr=tpr, fpr=fpr))
# }
