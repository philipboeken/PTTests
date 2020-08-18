experiment_sachs_lcd <- function(path = 'output/sachs/', observational = 1:8) {
  # Alternatively: observational = 1

  # Setup
  ##############################################

  context_vars <- c('AKT inh', 'G0076', 'Psitectorigenin', 'LY294002', 'U0126',
                    'PMA + noCD3/28', 'b2CAMP + noCD3/28')
  context_vars <- if (all(observational == 1)) context_vars else c('CD3/28', context_vars)

  system_vars <- setdiff(colnames(sachs_data), c(context_vars, 'experiment', 'CD3/28'))

  CX_combos <- as.matrix(expand.grid(C = context_vars, X = system_vars))
  lcd_triples <- as.matrix(expand.grid(C = context_vars, X = system_vars, Y = system_vars))
  lcd_triples <- lcd_triples[which(lcd_triples[, 2] != lcd_triples[, 3]),]

  get_results <- function(test, sachs_data, CX_combos, lcd_triples, observational) {
    `%dopar%` <- foreach::`%dopar%`
    CX_test_results <- foreach::foreach(i = 1:nrow(CX_combos), .combine = rbind) %dopar% {
      CX_combo <- CX_combos[i,]
      C <- CX_combo[1]
      X <- CX_combo[2]

      rows <- which(sachs_data$experiment %in% observational | sachs_data[[C]] == 1)
      data <- sachs_data[rows, CX_combo]
      CX <- test(data[[C]], data[[X]])
      return(data.frame(C = C, X = X, CX = CX))
    }

    result <- foreach::foreach(i = 1:nrow(lcd_triples), .combine = rbind) %dopar% {
      lcd_triple <- lcd_triples[i,]
      C <- lcd_triple[1]
      X <- lcd_triple[2]
      Y <- lcd_triple[3]

      rows <- which(sachs_data$experiment %in% observational | sachs_data[[C]] == 1)
      data <- sachs_data[rows, lcd_triple]
      CX <- CX_test_results[which(CX_test_results$C == C & CX_test_results$X == X), 'CX']
      XY <- test(data[[X]], data[[Y]])
      CY_X <- test(data[[C]], data[[Y]], data[[X]])
      return(data.frame(C = C, X = X, Y = Y, CX = CX, XY = XY, CY_X = CY_X))
    }
    return(result)
  }


  # Do test
  ##############################################

  # cores <- parallel::detectCores()
  # cl <- parallel::makeForkCluster(cores[1] - 1)
  # doParallel::registerDoParallel(cl)
  doParallel::registerDoParallel()

  results <- list(
    pcor = get_results(.ppcor_wrapper, sachs_data, CX_combos, lcd_triples, observational),
    polyatree = get_results(.pt_wrapper, sachs_data, CX_combos, lcd_triples, observational)
  )


  # Process results
  ##############################################

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  save(results, file = sprintf("%s%s.Rdata", path, timestamp))

  name <- 'sachs_pcor'
  .output_graph(results$pcor, path, name,
                alpha1 = list(strong = 0.0001, substantial = 0.005, weak = 0.05),
                alpha2 = list(strong = 0.05, substantial = 0.05, weak = 0.05))
  system(sprintf('dot -Tpdf %s%s.dot -o %s%s.pdf', path, name, path, name))

  name <- 'sachs_polyatree'
  .output_graph(results$polyatree, path, name,
                alpha1 = list(strong = 0.09, substantial = 0.2, weak = 0.5),
                alpha2 = list(strong = 0.91, substantial = 0.8, weak = 0.5))
  system(sprintf('dot -Tpdf %s%s.dot -o %s%s.pdf', path, name, path, name))
}
