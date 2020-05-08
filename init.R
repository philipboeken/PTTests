# source conditional independence tests
source('indepTests/splineGCM.R')
source('indepTests/bayesCITest.R')
source('indepTests/CCIT.R')

source('helpers.R')

suppressWarnings(library(ppcor))
suppressWarnings(library(RCIT))
suppressWarnings(library(reticulate))
use_python('/usr/local/bin/python3')
.ccit <- import('CCIT')
.ccit <- .ccit$CCIT$CCIT
