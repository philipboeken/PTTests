# Polya Tree independence tests (PTTests)

This R package accompanies the paper

'A Bayesian Nonparametric Conditional Two-sample Test with an Application to Local Causal Discovery'

## Installation

To install this R package, make sure you have R (version >= 3.5.0) and RStudio installed. Open ``PTTests.Rproj`` in RStudio and execute:

```R
install.packages("devtools")
devtools::install_local("./")
library(PTTests)
```

## Contents

This package implements the following Pólya tree-based hypothesis tests:

- ``pt_d_sample_ci_test`` (Anonymised)
- ``pt_d_sample_test`` (Holmes et al., 2015)
- ``pt_continuous_ci_test`` (Teymur and Filippi, 2019)
- ``pt_continuous_independence_test`` (Filippi and Holmes, 2017).

These tests may be called through the wrapper method ``pt_ci_test``.

The experiments from the paper can be reproduced by running the ``experiments.R`` script, or alternatively by executing the following functions:

- ``experiment_lcd_compare_tests_roc(simulation = 'paper')``
- ``experiment_lcd_roc_curves(simulation = 'paper')``
- ``experiment_sachs_lcd()``

The output is saved to the ``/output`` folder. You may provide a custom output folder using the ``path`` argument when invoking one of the above methods.

## Graphviz dependency

In the ``experiment_sachs_lcd()`` method, we make use of graphviz to generate .pdf files from .dot files. To ensure this works correctly, make sure you [have graphviz installed](http://www.graphviz.org/download/), and have the cli command ``dot`` working.

## Python dependency

In the ``experiment_lcd_compare_tests_roc()`` method we don't invoke the CCIT test by default. If you wish to obtain results for this test, some steps have to be taken. We run the CCIT test using the [provided python package](https://github.com/rajatsen91/CCIT), and approach this python package in R using ``reticulate``. For this to work, it is required to have Python 3 installed, and have the cli command ``python`` working. To install the CCIT package, run ``pip install CCIT == 0.4`` or ``sudo -H pip install CCIT == 0.4``.

Then uncomment from ``experiment_lcd_compare_tests``:

```R
  ccit = get_results(data, .ccit_wrapper),
```

and reinstall the package by running from an R console:

```R
devtools::install_local("./", force = TRUE)
```

## References

Filippi, S. and Holmes, C. C. (2017). A Bayesian nonparametric approach to testing for dependence between random variables. <em>Bayesian Analysis</em>, 12(4):919–938.

Holmes, C. C., Caron, F., Griffin, J. E., and Stephens, D. A. (2015). Two-sample Bayesian nonparametric hypothesis testing. <em>Bayesian Analysis</em>, 10(2):297–320.

Teymur, O. and Filippi, S. (2019). A Bayesian nonparametric test for conditional independence. <em>arXiv.org preprint</em>, arxiv:1910.11219 [stat.ME].
