
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intervalcalc

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/medewitt/intervalcalc/branch/main/graph/badge.svg)](https://codecov.io/gh/medewitt/intervalcalc?branch=main)
[![R-CMD-check](https://github.com/medewitt/intervalcalc/workflows/R-CMD-check/badge.svg)](https://github.com/medewitt/intervalcalc/actions)
<!-- badges: end -->

The goal of intervalcalc is to …

## Installation

You can install the released version of intervalcalc from
[r-universe](https://medewitt.r-universe.dev/ui#builds) with

``` r
options(repos = c(
    medewitt = 'https://medewitt.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
install.packages("intervalcalc")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("medewitt/intervalcalc")
```

You will also need to have a working installation of CmdStanR. This can
be downloaded and installed using:

``` r
install.packages("cmdstanr", 
    repos = c("https://mc-stan.org/r-packages/", 
    getOption("repos")))
```

Then you can install CmdStan using the following code from within a new
R session:

``` r
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan()
```

Full details are available at [the CmdStanR
website](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

## Code of Conduct

Please note that the intervalcalc project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
