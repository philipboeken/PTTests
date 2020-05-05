suppressMessages(library(ggplot2))
suppressMessages(library(plotly))

pplot <- function(x, y, count=FALSE) {
  if (count) {
    ggplot(data.frame(y,x), aes(x, y)) + geom_count()
  } else {
    result <- bayes.UCItest(x, y)
    cat("\nBF(H0, H1):\t", result$bf, "\nP(H1 | XY):\t", result$p_H1, "\n")
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() +
      labs(x = substitute(x), y = substitute(y))
  }
}

hhist <- function(x) {
  ggplot(data.frame(x), aes(x)) + geom_histogram(bins=length(x)/50)
}

gamplot <- function(x, y, binomial=FALSE) {
  if (binomial)
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x), method.args = list(family = "binomial"))
  else
    ggplot(data.frame(y,x), aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x))
}

pplot3d <- function(X,Y,Z) {
  scene <- list(
    xaxis = list(title = substitute(X)),
    yaxis = list(title = substitute(Y)),
    zaxis = list(title = substitute(Z))
  )
  fig <- plot_ly(x=X, y=Y, z=Z, type="scatter3d", mode="markers", color=Z)
  fig <- fig %>% layout(scene=scene)
  fig
}

pplotcs <- function(X, Y, Z, Z_min, Z_max) {
  XYZ <- cbind(X, Y, Z)
  data <- as.matrix(XYZ[Z_min < XYZ[,3] & XYZ[,3] < Z_max, c(1, 2)])
  pplot(data[,1], data[,2]) +
    labs(x = substitute(X), y = substitute(Y))
}

scale_norm <- function(X) {
  return(pnorm(X, mean=mean(X), sd=sd(X)))
}

scale_lin <- function(X) {
  if (ncol(as.matrix(X)) > 1) {
    X <- apply(X, 2, scale_lin)
  }
  return((X - min(X)) / (max(X) - min(X)))
}

