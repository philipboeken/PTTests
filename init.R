# source conditional independence tests
source('independence_tests/splineGCM.R')
source('independence_tests/bayesCITest.R')
source('independence_tests/CCIT.R')

source('helpers.R')

suppressWarnings(library(ppcor))
suppressWarnings(library(RCIT))
suppressWarnings(library(reticulate))
use_python('/usr/local/bin/python3')
.ccit <- import('CCIT')
.ccit <- .ccit$CCIT$CCIT
