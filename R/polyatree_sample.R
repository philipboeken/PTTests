polyatree_sample <- function(c = 1, depth = 12, aj = function(j) c*j^2, save = FALSE){
  
  ps <- list()
  ps[[1]] <- c(1)
  
  for(j in 1:depth) {
    alpha <- aj(j)
    theta <- DirichletReg::rdirichlet(2**(j-1), c(alpha, alpha))
    ind <- 2*seq_len(2**(j-1)) - 1
    ps[[j+1]] <- rep(0, 2**j)
    ps[[j+1]][ind] <- ps[[j]]
    ps[[j+1]][ind+1] <- ps[[j]]
    ps[[j+1]] <- ps[[j+1]] * as.vector(t(theta))
  }
  
  qs <- qnorm((1:(2**depth) - 0.5) / (2**depth))
  tmp <- qnorm((0:(2**depth)) / (2**depth))
  widths <- tmp[2:length(tmp)] - tmp[1:(length(tmp)-1)]
  y <- ps[[depth+1]] / widths
  
  plt <- ggplot2::ggplot(data=data.frame(x=qs, y=y), ggplot2::aes(x=x, y=y, group=1)) +
    ggplot2::geom_line() +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())
  if(save) {
    .ggsave('output/polyatree_sample', plt, 8, 8)
  }
  plt
}

