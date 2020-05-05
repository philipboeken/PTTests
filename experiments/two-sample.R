source('init.R', chdir=TRUE)
library(foreach)
library(doParallel)
library(distr)
library(ROCR)

gms <- function(n, C, theta) {
  # Gaussian mean shift
  (1-C) * rnorm(n) + C * rnorm(n, theta, 1)
}
gvs <- function(n, C, theta) {
  # Gaussian variance shift
  (1-C) * rnorm(n) + C * rnorm(n, 0, 1+theta)
}
gm <- function(n, C, theta) {
  # Gaussian mixture
  (1-C) * rnorm(n) + C * r(UnivarMixingDistribution(Norm(mean=-theta, sd=1), 
                                                    Norm(mean=theta, sd=1),
                                                    mixCoeff=c(0.5, 0.5)))(n)
}
t <- function(n, C, theta) {
  # Tails
  (1-C) * rnorm(n) + C * rt(n, 10^theta)
}
lms <- function(n, C, theta) {
  # Lognormal mean shift
  (1-C) * rlnorm(n) + C * rlnorm(n, theta, 1)
}
lvs <- function(n, C, theta) {
  # Lognormal variance shift
  (1-C) * rlnorm(n) + C * rlnorm(n, 0, 1+theta)
}

get_two_sample_data <- function(n) {
  theta <- sample(0:3, 1)
  p <- runif(1, 0.05, 0.95)
  C <- rbinom(n, 1, p)
  opts <- c(
    gms,
    gvs,
    gm,
    t
    # lms
    # lvs
  )
  X <- sample(opts, 1)[[1]](n, C, theta)
  
  data <- cbind(X, C)
  X1 <- data[data[,2] == 0, 1]
  X2 <- data[data[,2] == 1, 1]
  return(list(p=p, C=C, X=X, X1=X1, X2=X2, true=as.numeric(theta==0)))
  
}

get_results <- function(n, m){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_two_sample_data(n)
    return(data.frame(true=data$true,
             cor=cor.test(data$C, data$X)$p.value,
             # ks.test(data$X1, data$X2)$p.value,
             # kruskal.test(list(data$X1, data$X2))$p.value,
             # splineGCM(1, 2, c(), cbind(data$C, data$X)),
             RCoT=RCoT(data$C, data$X)$p,
             # RCIT(data$C, data$X)$p,
             BayesTS=bayes.UCItest(data$X, data$C, verbose=FALSE)$p_H0))
  }
  return(result)
}

cores <- detectCores()
cl <- makeForkCluster(cores[1]-1)
registerDoParallel(cl)
results <- get_results(500, 500)
stopCluster(cl)

roc_data <- c()
for (i in 2:ncol(results)) {
  pred <- prediction(results[,i], results[,1])
  res <- performance(pred, "tpr", "fpr")
  auc <- round(performance(pred, "auc")@y.values[[1]], 3)
  x <- res@x.values[[1]]
  y <- res@y.values[[1]]
  name <- colnames(results)[i]
  roc_data[[name]] <- list(data=data.frame(x=x, y=y), auc=auc, name=name)
}

plt <- ggplot()
for (roc in roc_data) {
  c <- paste(roc$name, ' (auc: ', roc$auc, ')', sep="")
  plt <- plt + geom_line(data=roc$data, aes(x, y, colour={{c}}))
}

plot(plt)


