CCIT <- function (X, Y, Z = NULL) {
  ## This is defined in test_wrappers.R to reduce overhead. Uncomment 
  ## when using this without test_wrappers.R.
  # library(reticulate)
  # use_python('/usr/local/bin/python3')
  # .ccit <- import('CCIT')
  # .ccit <- .ccit$CCIT$CCIT
  
  if (is.null(Z)){
    return(suppressWarnings(.ccit(matrix(X, ncol=1), matrix(Y, ncol=1), NULL)))
  }
  
  return(suppressWarnings(.ccit(matrix(X, ncol=1), matrix(Y, ncol=1), matrix(Z, ncol=1))))
}
