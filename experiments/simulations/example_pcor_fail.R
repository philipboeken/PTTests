source('independence_tests/test_wrappers.R')
source('experiments/simulations/maps.R')
source('helpers.R')
library(foreach)
library(doParallel)
suppressWarnings(library(cowplot))


# Input parameters
##############################################
n <- 300
m <- 500

err_sd <- 0.1

p_link <- 0.5
p_two_sample <- 0.5

nonlin_options <- c(parabolic)

interv_options <- c(mean_shift)



# Setup test
##############################################
get_data <- function(n, p_two_sample, p_link, err_sd, 
                     nonlin_options, interv_options, link=0) {
  C <- rbinom(n, 1, p_two_sample)
  
  # C -> Z <- X
  X <- rnorm(n)
  
  if (link < 0) {
    link_nonlin <- 0
  } else if (link == 0) {
    link_nonlin <- rbinom(1, 1, p_link)
  } else if (link > 0) {
    link_nonlin <- 1
  }
  
  Z <- link_nonlin * nonlin(nonlin_options, X)
  Z <- Z + err_sd * rnorm(n, 0, ifelse(sd(Z) > 0, sd(Z), 1/err_sd))
  
  Z <- do_intervention(interv_options, Z, C)
  
  cond_indep <- as.numeric(!link_nonlin)
  
  return(list(C=C, X=X, Z=Z, label=cond_indep))
}

get_results <- function(n, m, p_two_sample, p_link, err_sd,
                        nonlin_options, interv_options){
  result <- foreach(i=1:m, .combine=rbind) %dopar% {
    data <- get_data(n, p_two_sample, p_link, err_sd, nonlin_options, interv_options)
    return(data.frame(
      label=data$label,
      pcor=.pcor_wrapper(data$C, data$X, data$Z),
      bayes=.bayes_wrapper(data$C, data$X, data$Z)
    ))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- get_results(n, m, p_two_sample, p_link, err_sd,
                       nonlin_options, interv_options)

data_no_link <- get_data(n, p_two_sample, p_link, err_sd, nonlin_options, interv_options, -1)
data_linked <- get_data(n, p_two_sample, p_link, err_sd, nonlin_options, interv_options, 1)

stopCluster(.cl)


# Process results
##############################################
no_link <- data.frame(C=as.factor(data_no_link$C), X=data_no_link$X, Z=data_no_link$Z)
linked <- data.frame(C=as.factor(data_linked$C), X=data_linked$X, Z=data_linked$Z)

scat_plot_no_link <- ggplot() + 
  geom_point(data=no_link, aes(x=X, y=Z, colour=C)) + 
  theme(legend.title = element_blank(),
        legend.position = c(0.87, 0.12)) +
  scale_color_manual(labels = c("C=0", "C=1"), values = c("blue", "red"))

scat_plot_linked <- ggplot() + 
  geom_point(data=linked, aes(x=X, y=Z, colour=C)) + 
  theme(legend.title = element_blank(),
        legend.position = c(0.87, 0.12)) +
  scale_color_manual(labels = c("C=0", "C=1"), values = c("blue", "red"))

roc_plot <- pplot_roc(results[,1], results[,-1], NULL, c(0.8, 0.12))

grid <- plot_grid(scat_plot_no_link, scat_plot_linked, roc_plot, nrow=1)
plot(grid)

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
.path <- 'experiments/simulations/output/example-pcor-fail/'

save.image(file=paste(.path, timestamp, '.Rdata', sep=''))

.ggsave(paste(.path, timestamp, sep=''), grid, 30, 10)
.ggsave(paste(.path, 'last', sep=''), grid, 30, 10)
