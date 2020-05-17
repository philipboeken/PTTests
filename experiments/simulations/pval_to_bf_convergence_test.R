# Clear workspace
# rm(list = ls(all.names = TRUE))
# gc()

# Imports
source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('independence_tests/RhoBFP.R')
source('helpers.R')
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))
suppressWarnings(library(cowplot))
require(scales)


.data_1 <- function(n) {
  return(list(X=runif(n), Y=runif(n), Z=runif(n)))
}

.data_2 <- function(n) {
  U <- runif(n)
  V <- runif(n)
  Z <- runif(n)
  return(list(X=(U+Z)/2, Y=(V+Z)/2, Z=Z))
}

.data_3 <- function(n) {
  X <- runif(n)
  Y <- runif(n)
  W <- rnorm(n, 0, 1/2)
  Z <- atan(Y / X + W) / pi + 1/2
  return(list(X=X, Y=Y, Z=Z))
}

.data_4 <- function(n) {
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

get_results <- function(nr, n, m, test){
  get_data <- switch(nr, .data_1, .data_2, .data_3, .data_4)
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n)
    return(c(
      test(data$X, data$Y) 
      # 1-test(data$X, data$Y, data$Z)
    ))
  }
  return(result)
}

Ns <- c(3, 5, 10, 21, 46, 100, 215, 464, 1000)
m <- 150

.do_tests <- function(test) {
  results <- list()
  for (nr in 1:4) {
    var <- paste('test', nr, sep="")
    results[[var]] <- list(uci=c(), ci=c())
    for (n in Ns) {
      cat(format(Sys.time(), "%X"), "Running test", nr, "with n =", n, "\n")
      res <- get_results(nr, n, m, test)
      results[[var]]$uci <- rbind(results[[var]]$uci, quantile(res, c(0.05, 0.25, 0.5, 0.75, 0.95)))
      # results[[var]]$ci <- rbind(results[[var]]$ci, quantile(res[,2], c(0.05, 0.25, 0.5, 0.75, 0.95)))
    }
  }
  return(results)
}

# Setup parallel backend to use many processors
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

# Perform tests:

results <- list(
  pcor=.do_tests(.pcor_wrapper),
  bayes=.do_tests(.bayes_wrapper),
  # bcor_wg=.do_tests(.bcor_wg_wrapper),
  pcor_bayes_ly=.do_tests(.bayes_transform2(.pcor_wrapper)),
  pcor_bayes_pb=.do_tests(.bayes_transform(.pcor_wrapper)),
  rcot=.do_tests(.rcot_wrapper),
  rcot_bayes_pb=.do_tests(.bayes_transform(.rcot_wrapper))
  # rcot_bayes_ly=.do_tests(.bayes_transform2(.rcot_wrapper))
  # gcm_freq2bayes=.do_tests(.bayes_transform(.gcm_wrapper))
)


stopCluster(.cl)

# Generate plots:
.get_plot <- function(data, title) {
  return(ggplot(data.frame(cbind(data, Ns)), aes(Ns)) +
           geom_line(aes(y=X5., colour="5")) +
           geom_line(aes(y=X25., colour="25")) +
           geom_line(aes(y=X50.)) +
           geom_line(aes(y=X75., colour="25")) +
           geom_line(aes(y=X95., colour="5")) +
           scale_x_continuous(limits=c(1, max(Ns)), trans = log10_trans()) +
           ylim(0, 1) +
           labs(x = "N", y="P(X indep Y | data)", title=title) +
           theme(legend.position = "none",
                 plot.title = element_text(size=12, hjust=0.5))
  )
}

grid <- plot_grid(
  .get_plot(results$pcor$test1$uci, "pcor (X indep Y)"),
  .get_plot(results$pcor$test2$uci, "pcor (X dep Y)"),
  # .get_plot(results$pcor$test3$uci, "pcor Test 3 UCI"),
  # .get_plot(results$pcor$test4$uci, "pcor Test 4 UCI"),
  .get_plot(results$bayes$test1$uci, "Bayes (X indep Y)"),
  .get_plot(results$bayes$test2$uci, "Bayes (X dep Y)"),
  # .get_plot(results$bayes$test3$uci, "Bayes Test 3 UCI"),
  # .get_plot(results$bayes$test4$uci, "Bayes Test 4 UCI"),
  # .get_plot(results$bcor_wg$test1$uci, "bcor_wg (X indep Y)"),
  # .get_plot(results$bcor_wg$test2$uci, "bcor_wg (X dep Y)"),
  # .get_plot(results$bcor_wg$test3$uci, "bcor_wg Test 3 UCI"),
  # .get_plot(results$bcor_wg$test4$uci, "bcor_wg Test 4 UCI"),
  .get_plot(results$pcor_bayes_ly$test1$uci, "bcor_ly (X indep Y)"),
  .get_plot(results$pcor_bayes_ly$test2$uci, "bcor_ly (X dep Y)"),
  # .get_plot(results$pcor_bayes_ly$test3$uci, "bcor_ly Test 3 UCI"),
  # .get_plot(results$pcor_bayes_ly$test4$uci, "bcor_ly Test 4 UCI"),
  .get_plot(results$pcor_bayes_pb$test1$uci, "bcor_pb (X indep Y)"),
  .get_plot(results$pcor_bayes_pb$test2$uci, "bcor_pb (X dep Y)"),
  # .get_plot(results$pcor_bayes_pb$test3$uci, "bcor_pb Test 3 UCI"),
  # .get_plot(results$pcor_bayes_pb$test4$uci, "bcor_pb Test 4 UCI"),
  .get_plot(results$rcot$test1$uci, "rcot (X indep Y)"),
  .get_plot(results$rcot$test2$uci, "rcot (X dep Y)"),
  # .get_plot(results$rcot$test3$uci, "rcot Test 3 UCI"),
  # .get_plot(results$rcot$test4$uci, "rcot Test 4 UCI"),
  .get_plot(results$rcot_bayes_pb$test1$uci, "rcot_bayes (X indep Y)"),
  .get_plot(results$rcot_bayes_pb$test2$uci, "rcot_bayes (X dep Y)"),
  # .get_plot(results$rcot_bayes_pb$test3$uci, "rcot_bayes Test 3 UCI"),
  # .get_plot(results$rcot_bayes_pb$test4$uci, "rcot_bayes Test 4 UCI"),
  # .get_plot(results$rcot_bayes_ly$test1$uci, "rcot_bayes_ly (X indep Y)"),
  # .get_plot(results$rcot_bayes_ly$test2$uci, "rcot_bayes_ly (X dep Y)"),
  # .get_plot(results$rcot_bayes_ly$test3$uci, "rcot_bayes_ly Test 3 UCI"),
  # .get_plot(results$rcot_bayes_ly$test4$uci, "rcot_bayes_ly Test 4 UCI"),
  nrow=6
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

.path <- 'experiments/simulations/output/pval-to-bf-convergence-tests/'

save.image(file=paste(.path, timestamp, ".Rdata", sep=""))

.ggsave(paste(.path, timestamp, sep=""), grid, 30, 35)
.ggsave(paste(.path, 'last', sep=""), grid, 20, 35)

# plot(grid)
