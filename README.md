<!-- badges: starts -->
[![R-CMD-check](https://github.com/resplab/predtools/workflows/R-CMD-check/badge.svg)](https://github.com/resplab/predtools/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/predtools)](https://cran.r-project.org/package=predtools)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/predtools)](https://cran.r-project.org/package=predtools)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

# Overview

`predtools` provides miscellaneous tools for developing and evaluating
prediction models.

# Table of Contents

-   [Installation](#installation)
-   [Example](#example)
-   [Model-based
    ROC](https://resplab.github.io/predtools/articles/mROC.html)
-   [Intercept
    Adjustment](https://resplab.github.io/predtools/articles/interceptAdj.html)
-   [Calibration
    Plot](https://resplab.github.io/predtools/articles/calibPlot.html)
-   [Unit Normal Loss Integral in Two Dimensions](https://resplab.github.io/predtools/articles/UNLI2D.html)

## Installation

You can install the released version of predtools from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("predtools")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("resplab/predtools")
```

## Example

The function `calibration_plot` takes observed and predicted values from
a prediction model and uses ggplot2 to produce a calibration plot:

    library(predtools)
    library(dplyr)
    x <- rnorm(100, 10, 2)
    y <- x + rnorm(100,0, 1)
    data <- tibble(x,y)
    calibration_plot(data, obs = "x", pred_1 = "y")

See vignettes for more advanced functionalities, including [model-based
ROC](https://resplab.github.io/predtools/articles/mROC.html), [intercept
adjustment](https://resplab.github.io/predtools/articles/interceptAdj.html), [calibration
plot](https://resplab.github.io/predtools/articles/calibPlot.html), and [unit normal loss integral in two dimensions](https://resplab.github.io/predtools/articles/UNLI2D.html).

You can also access the vignettes from R:

    browseVignettes("predtools")
