experiment_teymur_filippi <- function(path = 'output/teymur_filippi/') {
  library(cowplot)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  require(scales)
  
  data_1 <- function(n) {
    return(list(X=runif(n), Y=runif(n), Z=runif(n)))
  }
  
  data_2 <- function(n) {
    U <- runif(n)
    V <- runif(n)
    Z <- runif(n)
    return(list(X=(U+Z)/2, Y=(V+Z)/2, Z=Z))
  }
  
  data_3 <- function(n) {
    X <- runif(n)
    Y <- runif(n)
    W <- rnorm(n, 0, 1/2)
    Z <- atan(Y / X + W) / pi + 1/2
    return(list(X=X, Y=Y, Z=Z))
  }
  
  data_4 <- function(n) {
    suppressWarnings(library(mvtnorm))
    mv <- rmvnorm(n, mean=c(0.5, 0.5), sigma=cbind(c(2, 1), c(1, 2))/20)
    mv[((mv[,1] < 0) | (1 < mv[,1])) & ((mv[,2] < 0) | (1 < mv[,2])), ] <- 0
    X <- mv[,1]
    Y <- mv[,2]
    B <- rbinom(n, 1, 0.9)
    Z <- rnorm(n, 1/2, 1/2)
    Z[(Z < 0) | (1 < Z)] <- 0
    Z <- Z * (B + (1-B) * X * Y)
    return(list(X=X, Y=Y, Z=Z))
  }
  
  get_results <- function(nr, n){
    get_data <- switch(nr, data_1, data_2, data_3, data_4)
    result <- foreach(i=1:100, .combine=rbind) %dopar% {
      data <- get_data(n)
      return(c(polyatree_continuous_independence_test(data$X, data$Y)$p_H1,
               polyatree_continuous_ci_test(data$X, data$Y, data$Z, rho=0)$p_H1))
    }
    return(result)
  }
  
  Ns <- c(1, 5, 10, seq(20, 100, by=20), seq(200, 500, by=100), 600, 800, 1000)
  # Ns <- c(50, 100)
  
  do_test <- function(nr) {
    result <- list(uci=c(), ci=c())
    for (n in Ns) {
      cat(format(Sys.time(), "%X"), "Running test", nr, "with n =", n, "\n")
      res <- get_results(nr, n)
      result$uci <- rbind(result$uci, quantile(res[,1], c(0.05, 0.25, 0.5, 0.75, 0.95)))
      result$ci <- rbind(result$ci, quantile(res[,2], c(0.05, 0.25, 0.5, 0.75, 0.95)))
    }
    return(result)
  }
  
  # Setup parallel backend to use many processors
  cores <- detectCores()
  cl <- makeForkCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Perform tests:
  t1 <- do_test(1)
  t2 <- do_test(2)
  t3 <- do_test(3)
  t4 <- do_test(4)
  
  stopCluster(cl)
  
  # Generate plots:
  get_plot <- function(data, title) {
    return(ggplot(data.frame(cbind(data, Ns)), aes(Ns)) +
             geom_line(aes(y=X5., colour="5")) +
             geom_line(aes(y=X25., colour="25")) +
             geom_line(aes(y=X50.)) +
             geom_line(aes(y=X75., colour="25")) +
             geom_line(aes(y=X95., colour="5")) +
             scale_x_continuous(limits=c(1, max(Ns)), trans = log10_trans()) +
             ylim(0, 1) +
             labs(x = "N", y="P(H1 | W)", title=title) +
             theme(legend.position = "none",
                   plot.title = element_text(size=12, hjust=0.5))
    )
  }
  
  grid <- plot_grid(get_plot(t1$uci, "Test 1, MI"),
                    get_plot(t2$uci, "Test 2, MI"),
                    get_plot(t3$uci, "Test 3, MI"),
                    get_plot(t4$uci, "Test 4, MI"),
                    get_plot(t1$ci, "Test 1, CI"),
                    get_plot(t2$ci, "Test 2, CI"),
                    get_plot(t3$ci, "Test 3, CI"),
                    get_plot(t4$ci, "Test 4, CI"),
                    nrow=2)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  save.image(file=paste(path, 'teymur_filippi_', timestamp, ".Rdata", sep=""))
  
  ggsave(
    paste(path, 'teymur_filippi_', timestamp, ".pdf", sep=""),
    plot = grid,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 30,
    height = 15,
    units = "cm",
    dpi = 300,
    limitsize = TRUE
  )
}
