source('independence_tests/test_wrappers.R')
source('helpers.R')
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))
suppressWarnings(library(cowplot))
suppressWarnings(library(readr))


# Input parameters
##############################################
test_ensembles <- c(
  .bayes_wrapper,
  .pcor_wrapper,
  .gcm_wrapper,
  .rcot_wrapper,
  .ccit_wrapper
)


# Setup test
##############################################
.col_types <- cols('AKT inh'=col_integer(), 'G0076'=col_integer(),
                   'LY294002'=col_integer(), 'PMA/beta2CAMP + noAlphaCD3/28'=col_integer(), 
                   'Psitectorigenin'=col_integer(), 'U0126'=col_integer(),
                   'experiment'=col_integer())

sachs_data_pooled <- read_csv("experiments/sachs/sachs_data_pooled.csv", 
                              col_types = .col_types)

.context_vars <- c('AKT inh', 'G0076', 'LY294002', 'PMA/beta2CAMP + noAlphaCD3/28', 
                      'Psitectorigenin', 'U0126', 'experiment')

.system_vars <- setdiff(colnames(sachs_data_pooled), .context_vars)

pos_lcd_triples <- as.matrix(expand.grid(C=.context_vars, X=.system_vars, Y=.system_vars))
pos_lcd_triples <- pos_lcd_triples[which(pos_lcd_triples[ ,2] != pos_lcd_triples[ ,3]), ]

get_results <- function(test) {
  result <- foreach(i=1:nrow(pos_lcd_triples), .combine=rbind) %dopar% {
    lcd_triple <- pos_lcd_triples[i,]
    data <- as.matrix(sachs_data_pooled[, lcd_triple])
    ts <- test(data[,1], data[,2])
    uci <- test(data[,2], data[,3])
    ci <- test(data[,1], data[,2], data[,3])
    lcd <- (1 - ts) * (1 - uci) * ci
    return(data.frame(ts=ts, uci=uci, ci=ci, lcd=lcd))
  }
  return(result)
}


# Do test
##############################################
.cores <- detectCores()
.cl <- makeForkCluster(.cores[1]-1)
registerDoParallel(.cl)

results <- list(
  pcor=get_results(.pcor_wrapper),
  bayes=get_results(.bayes_wrapper),
  gcm=get_results(.gcm_wrapper),
  ccit=get_results(.ccit_wrapper),
  rcot=get_results(.rcot_wrapper)
)

stopCluster(.cl)


# Process results
##############################################

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save.image(file=paste('experiments/sachs/output/sachs_output_', timestamp, ".Rdata", sep=""))
