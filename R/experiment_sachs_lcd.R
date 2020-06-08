experiment_sachs_lcd <- function (path = 'R/output/sachs/') {
  library(doParallel)
  
  # Setup test
  ##############################################
  col_types <- readr::cols('AKT inh' = readr::col_integer(), 'G0076' = readr::col_integer(),
                           'LY294002' = readr::col_integer(), 'PMA + noCD3/28' = readr::col_integer(), 
                           'b2CAMP + noCD3/28' = readr::col_integer(), 
                           'Psitectorigenin' = readr::col_integer(), 
                           'U0126' = readr::col_integer(), 
                           # 'CD3/28' = readr::col_integer(), 
                           'experiment' = readr::col_integer())
  
  sachs_data_pooled <- readr::read_csv("R/data/sachs_data.csv", col_types = col_types)
  
  context_vars <- c('AKT inh', 'G0076', 'LY294002', 'PMA + noCD3/28', 'b2CAMP + noCD3/28', 
                    # 'CD3/28',
                    'Psitectorigenin', 'U0126')
  
  system_vars <- setdiff(colnames(sachs_data_pooled), c(context_vars, 'experiment'))
  
  CX_combos <- as.matrix(expand.grid(C = context_vars, X = system_vars))
  lcd_triples <- as.matrix(expand.grid(C = context_vars, X = system_vars, Y = system_vars))
  lcd_triples <- lcd_triples[which(lcd_triples[ ,2] != lcd_triples[ ,3]), ]
  
  get_results <- function(test, sachs_data_pooled, CX_combos, lcd_triples) {
    CX_test_results <- foreach::foreach(i = 1:nrow(CX_combos), .combine = rbind) %dopar% {
      CX_combo <- CX_combos[i,]
      C <- CX_combo[1]
      X <- CX_combo[2]
      
      data <- sachs_data_pooled[, CX_combo]
      CX <- test(data[[C]], data[[X]])
      return(data.frame(C = C, X = X, CX = CX))
    }
    
    result <- foreach::foreach(i = 1:nrow(lcd_triples), .combine = rbind) %dopar% {
      lcd_triple <- lcd_triples[i,]
      C <- lcd_triple[1]
      X <- lcd_triple[2]
      Y <- lcd_triple[3]
      
      data <- sachs_data_pooled[, lcd_triple]
      CX <- CX_test_results[which(CX_test_results$C == C &  CX_test_results$X == X), 'CX']
      XY <- test(data[[X]], data[[Y]])
      CY_X <- test(data[[C]], data[[Y]], data[[X]])
      return(data.frame(C = C, X = X, Y = Y, CX = CX, XY = XY, CY_X = CY_X))
    }
    return(result)
  }
  
  
  # Do test
  ##############################################
  cores <- detectCores()
  cl <- makeForkCluster(cores[1]-1)
  registerDoParallel(cl)
  
  results <- list(
    pcor = get_results(.ppcor_wrapper, sachs_data_pooled, CX_combos, lcd_triples),
    polyatree = get_results(.polyatree_wrapper, sachs_data_pooled, CX_combos, lcd_triples)
  )
  
  stopCluster(cl)
  
  
  # Process results
  ##############################################
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  name <- 'pcor'
  .output_graph(results$pcor, path, name, 
                alpha1 = list(strong = 0.0001, substantial = 0.005, weak = 0.05),
                alpha2 = list(strong = 0.05, substantial = 0.05, weak = 0.05))
  system(sprintf('dot -Tpdf %s%s.dot -o %s%s.pdf', path, name, path, name))
  
  name <- 'polyatree'
  .output_graph(results$polyatree, path, name, 
                alpha1 = list(strong = 0.09, substantial = 0.2, weak = 0.5),
                alpha2 = list(strong = 0.91, substantial = 0.8, weak = 0.5))
  system(sprintf('dot -Tpdf %s%s.dot -o %s%s.pdf', path, name, path, name))
}
