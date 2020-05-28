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

obs <- FALSE
no_CD <- FALSE

# Setup test
##############################################
.col_types <- cols('AKT inh'=col_integer(), 'G0076'=col_integer(),
                   'LY294002'=col_integer(), 'PMA'=col_integer(), 'b2CAMP'=col_integer(), 
                   'Psitectorigenin'=col_integer(), 'U0126'=col_integer(),
                   'experiment'=col_integer())

sachs_data_pooled <- read_csv("experiments/sachs/sachs_data.csv", 
                              col_types = .col_types)

.context_vars <- c('AKT inh', 'G0076', 'LY294002', 'PMA', 'b2CAMP', 'Psitectorigenin', 'U0126')

if (no_CD) {
  drop <- c('PMA', 'b2CAMP')
  sachs_data_pooled <- sachs_data_pooled[, !(names(sachs_data_pooled) %in% drop)]
  .context_vars <- c('AKT inh', 'G0076', 'LY294002', 'Psitectorigenin', 'U0126')
}

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
    CX_combo <- CX_combos[i,]
    C <- CX_combo[1]
    X <- CX_combo[2]
    
    idx <- TRUE
    if (obs) {
      idx <- which(sachs_data_pooled[[C]] == 1 | sachs_data_pooled$experiment == 1)
    }
    
    data <- sachs_data_pooled[idx, CX_combo]
    CX <- test(data[[C]], data[[X]])
    return(data.frame(C=C, X=X, CX=CX))
  }
  
  result <- foreach(i=1:nrow(lcd_triples), .combine=rbind) %dopar% {
    lcd_triple <- lcd_triples[i,]
    C <- lcd_triple[1]
    X <- lcd_triple[2]
    Y <- lcd_triple[3]
    
    idx <- TRUE
    if (obs) {
      idx <- which(sachs_data_pooled[[C]] == 1 | sachs_data_pooled$experiment == 1)
    }
    
    data <- sachs_data_pooled[idx, lcd_triple]
    CX <- CX_test_results[which(CX_test_results$C == C &  CX_test_results$X == X), 'CX']
    XY <- test(data[[X]], data[[Y]])
    CY_X <- test(data[[C]], data[[Y]], data[[X]])
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

name <- 'bayes_weak'
if (obs) {
  name <- paste(name, '_obs', sep="")
}
if (no_CD) {
  name <- paste(name, '_no_CD', sep="")
}
bayes_triples_weak <- filter(results$bayes, CX <= 0.5, XY <= 0.5, CY_X >= 0.5)
.output_graph(bayes_triples_weak, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_weak, paste(.path, name, sep=''))

name <- 'bayes_substantial'
if (obs) {
  name <- paste(name, '_obs', sep="")
}
if (no_CD) {
  name <- paste(name, '_no_CD', sep="")
}
bayes_triples_substantial <- filter(results$bayes, CX <= 0.2, XY <= 0.2, CY_X >= 0.8)
.output_graph(bayes_triples_substantial, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_substantial, paste(.path, name, sep=''))

name <- 'bayes_strong'
if (obs) {
  name <- paste(name, '_obs', sep="")
}
if (no_CD) {
  name <- paste(name, '_no_CD', sep="")
}
bayes_triples_strong <- filter(results$bayes, CX <= 0.09, XY <= 0.09, CY_X >= 0.91)
.output_graph(bayes_triples_strong, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
.output_sachs_plots(sachs_data_pooled, bayes_triples_strong, paste(.path, name, sep=''))

name <- 'bayes_levels'
if (obs) {
  name <- paste(name, '_obs', sep="")
}
if (no_CD) {
  name <- paste(name, '_no_CD', sep="")
}
.output_graph_levels(bayes_triples_strong, bayes_triples_substantial, bayes_triples_weak, .path, name)
system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))

# name <- 'pcor'
# if (obs) {
#   name <- paste(name, '_obs', sep="")
# }
# if (no_CD) {
#   name <- paste(name, '_no_CD', sep="")
# }
# pcor_triples <- filter(results$pcor, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(pcor_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, pcor_triples, paste(.path, name, sep=''))
# 
# name <- 'rcot'
# if (obs) {
#   name <- paste(name, '_obs', sep="")
# }
# if (no_CD) {
#   name <- paste(name, '_no_CD', sep="")
# }
# rcot_triples <- filter(results$rcot, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(rcot_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, rcot_triples, paste(.path, name, sep=''))
# 
# name <- 'gcm'
# if (obs) {
#   name <- paste(name, '_obs', sep="")
# }
# if (no_CD) {
#   name <- paste(name, '_no_CD', sep="")
# }
# gcm_triples <- filter(results$gcm, CX <= 0.01, XY <= 0.01, CY_X >= 0.01)
# .output_graph(gcm_triples, .path, name)
# system(sprintf('dot -Tpdf %s%s.dot -o %s%s_graph.pdf', .path, name, .path, name))
# .output_sachs_plots(sachs_data_pooled, gcm_triples, paste(.path, name, sep=''))

save.image(file=paste(.path, timestamp, ".Rdata", sep=""))


