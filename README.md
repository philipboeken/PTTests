# Polya Tree independence tests

This R package accompanies the paper 'A Bayesian Nonparametric Conditional Two-sample Test with an Application to Local Causal Discovery'.

## Installation

To install this R package, make sure you have R installed (version >= 3.5.0), open an R console and execute

```R
install.packages("devtools")
library(devtools)
devtools::install_local("<full path to package>/polyatreelcd")
```

## Contents

This package implements the following Pólya tree based hypothesis tests:

- ``polyatree_two_sample_ci_test`` (this work)
- ``polyatree_two_sample_test`` (Holmes et al., 2015)
- ``polyatree_continuous_ci_test`` (Teymur and Filippi, 2019)
- ``polyatree_continuous_independence_test`` (Filippi and Holmes, 2017).

These tests may be called through the wrapper method ``polyatree_ci_test``.

The experiments from the paper can be reproduced by executing the following methods:

- ``experiment_pcor_fail``
- ``experiment_lcd_roc_curves``
- ``experiment_lcd_compare_tests``
- ``experiment_sachs_lcd``

The output is saved to the ``/output`` folder. You may provide a path to your custom output folder using the ``path`` argument when invoking one of the above methods.

## Graphviz dependency

We make use of graphviz to generate .pdf files from .dot files. To ensure this works correctly, make sure you [have graphviz installed](http://www.graphviz.org/download/), and have ``dot`` cli command working.

## Python dependency

By default, we don't invoke the CCIT test in the ``experiment_lcd_compare_tests`` method. If you wish to obtain results for this test, some steps have to be taken. We run the CCIT test using the [provided python package](https://github.com/rajatsen91/CCIT), and approach this python package in R using ``reticulate``. For this to work, it is required to have python3 installed, and have a binary located at ``/usr/local/bin/python3``. To install the CCIT package, run ``pip install CCIT == 0.4`` or ``sudo -H pip install CCIT == 0.4``.

Then uncomment from ``experiment_lcd_compare_tests``:

```R
  ccit = get_results(data, .ccit_wrapper),
```

and reinstall the package by running from an R console:

```R
devtools::install_local("<path to package>/polyatreelcd", force = TRUE)
```
