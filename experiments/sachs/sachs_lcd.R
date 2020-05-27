# Clear workspace
rm(list = ls(all.names = TRUE))

# Imports
source('independence_tests/test_wrappers.R')
source('helpers.R')
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))
suppressWarnings(library(cowplot))
suppressWarnings(library(readr))
suppressWarnings(library(Rgraphviz))


# Setup test
##############################################
.col_types <- cols('AKT inh'=col_integer(), 'G0076'=col_integer(),
                   'LY294002'=col_integer(), 'PMA/beta2CAMP + noAlphaCD3/28'=col_integer(), 
                   'Psitectorigenin'=col_integer(), 'U0126'=col_integer(),
                   'experiment'=col_integer())

sachs_data_pooled <- read_csv("experiments/sachs/sachs_data_pooled.csv", 
                              col_types = .col_types)

.context_vars <- c('AKT inh', 'G0076', 'LY294002', 'PMA/beta2CAMP + noAlphaCD3/28', 
                   'Psitectorigenin', 'U0126')

.system_vars <- setdiff(colnames(sachs_data_pooled), c(.context_vars, 'experiment'))

# Remove zero's?
# idx <- apply(sachs_data_pooled, 1, function(row) all(row[.system_vars] !=0 ))
# sachs_data_pooled <- sachs_data_pooled[idx,]

CX_combos <- as.matrix(expand.grid(C=.context_vars, X=.system_vars))
XY_combos <- as.matrix(expand.grid(X=.system_vars, Y=.system_vars))
XY_combos <- CX_combos[which(CX_combos[ ,1] != CX_combos[ ,2]), ]
lcd_triples <- as.matrix(expand.grid(C=.context_vars, X=.system_vars, Y=.system_vars))
lcd_triples <- lcd_triples[which(lcd_triples[ ,2] != lcd_triples[ ,3]), ]

get_results <- function(test) {
  CX_test_results <- foreach(i=1:nrow(CX_combos), .combine=rbind) %dopar% {
    C <- CX_combos[i,1]
    X <- CX_combos[i,2]
    idx <- which(sachs_data_pooled[, C] == 1 | sachs_data_pooled[,'experiment'] == 1)
    data <- as.matrix(sachs_data_pooled[idx, c(C, X)])
    CX <- test(data[,1], data[,2])
    return(data.frame(C=C, X=X, CX=CX))
  }
  
  result <- foreach(i=1:nrow(lcd_triples), .combine=rbind) %dopar% {
    lcd_triple <- lcd_triples[i,]
    C <- lcd_triple[1]
    X <- lcd_triple[2]
    Y <- lcd_triple[3]
    idx <- which(sachs_data_pooled[,C] == 1 | sachs_data_pooled[,'experiment'] == 1)
    data <- as.matrix(sachs_data_pooled[idx, lcd_triple])
    CX <- CX_test_results[which(CX_test_results[,'C'] == C &  CX_test_results[,'X'] == X), 'CX']
    XY <- test(data[,2], data[,3])
    CY_X <- test(data[,1], data[,3], data[,2])
    return(data.frame(C=C, X=X, Y=Y, CX=CX, XY=XY, CY_X=CY_X))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- list(
  # pcor=get_results(.pcor_wrapper),
  # rcot=get_results(.rcot_wrapper),
  # gcm=get_results(.gcm_wrapper),
  bayes=get_results(.bayes_wrapper)
)

stopCluster(.cl)


# Process results
##############################################

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
.path <- 'experiments/sachs/output/'

name <- 'sachs_output_bayes_weak'
bayes_triples_weak <- filter(results$bayes, CX <= 0.5, XY <= 0.5, CY_X >= 0.5)
.output_graph(bayes_triples_weak, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_weak, paste(.path, name, sep=''))

name <- 'sachs_output_bayes_substantial'
bayes_triples_substantial <- filter(results$bayes, CX <= 0.2, XY <= 0.2, CY_X >= 0.8)
.output_graph(bayes_triples_substantial, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_substantial, paste(.path, name, sep=''))

name <- 'sachs_output_bayes_strong'
bayes_triples_strong <- filter(results$bayes, CX <= 0.09, XY <= 0.09, CY_X >= 0.91)
.output_graph(bayes_triples_strong, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_strong, paste(.path, name, sep=''))

name <- 'sachs_output_bayes_levels'
.output_graph_levels(bayes_triples_strong, bayes_triples_substantial, bayes_triples_weak, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))

# name <- 'sachs_output_pcor'
# pcor_triples <- filter(results$pcor, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(pcor_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, pcor_triples, paste(.path, name, sep=''))
# 
# name <- 'sachs_output_rcot'
# rcot_triples <- filter(results$rcot, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(rcot_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, rcot_triples, paste(.path, name, sep=''))
# 
# name <- 'sachs_output_gcm'
# gcm_triples <- filter(results$gcm, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(gcm_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, gcm_triples, paste(.path, name, sep=''))

save.image(file=paste(.path, timestamp, ".Rdata", sep=""))


