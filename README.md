# Polya Tree independence tests

This R package accompanies the paper 'A Bayesian Nonparametric Conditional Two-sample Test with an Application to Local Causal Discovery'.

## Contents

The main file is ``independence_tests/polyatreeTests.R``, which implements the following:

- conditional two-sample test (this work)
- two-sample test (Holmes 2015)
- continuous conditional independence test (Teymur and Filippi, 2019)
- continuous independence test (Filippi and Holmes, 2017).

These tests should be called through the wrapper method ``polyatree.CITest``.

The experiments from the paper can be reproduced by sourcing from the project root:

- ``experiments/simulations/example_pcor_fail.R``
- ``experiments/simulations/lcd_triple_roc_test.R``
- ``experiments/simulations/lcd_measures_roc_test.R``
- ``experiments/sachs/sachs_lcd.R``

## Dependencies

This R package requires installation of the following R packages

- cowplot
- foreach
- doParallel
- GeneralisedCovarianceMeasure
- invgamma
- latex2exp
- matrixStats
- ppcor
- RCIT (via GitHub)
- readr
- reticulate

To install the RCIT package, run:

```R
library(devtools)
install_github("ericstrobl/RCIT")
```
