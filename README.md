Tutorial for model-based ROC (mROC) analysis
================

# predtools

<!-- badges: start -->

[![R-CMD-check](https://github.com/resplab/predtools/workflows/R-CMD-check/badge.svg)](https://github.com/resplab/predtools/actions)
<!-- badges: end -->

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

## What is predtools?

Imagine you have developed a risk prediction model using some
development dataset. The risk prediction model takes in some predictors
(e.g., sex, age, previous disease history) and returns the risk of an
event (e.g., risk of disease relapse in the next 12 months). You would
like to evaluate the performance of the risk model in a new (external)
validation sample. Among other things, you typically evaluate the
Receiver Operating Characteristic (ROC) curve of the risk prediciton
model in the new sample.

Now, model-based ROC (mROC) curve is the ROC curve that should be
observed if the prediction model is calibrated in the external
population. Comparing the empirical ROC and mROC curves in the new
sample can be informative on if the model is calibrated in the new
sample.

## How the package works

The package provides simple functions for mROC-related methods. It also
comes with exemplary datasets. Below we provide a step-by-step
illustration

## A step-by-step guide

Imagine the variable y indicates risk of disease recurrence in a unit of
time. We have a prediction model that quantifies this risk given a
patient’s age, disease severity level, sex, and whether the patient has
comorbidity.

The package comes with two exemplary datasets. dev\_data and val\_data.
We use the dev\_data as the development sample and the val\_data as the
external validation sample.

``` r
library(predtools)
data(dev_data)
data(val_data)
```

dev\_data has 500 rows. val\_data has 400 rows.

Here are the first few rows of dev\_data:

|        age |   severity | sex | comorbidity |   y |
|-----------:|-----------:|----:|------------:|----:|
|  0.5160477 |  0.2593806 |   0 |           1 |   1 |
|  0.1572608 | -1.2302576 |   0 |           0 |   0 |
|  1.2804923 |  1.8332095 |   1 |           1 |   1 |
|  1.1347625 |  1.3272407 |   0 |           0 |   1 |
|  0.8355657 |  0.3424236 |   0 |           0 |   1 |
|  0.3522012 |  2.1135965 |   0 |           1 |   1 |
| -0.5431004 | -2.9125308 |   0 |           1 |   0 |

We use the development data to fit a logistic regression model as our
risk prediction model:

``` r
reg<-glm(y~sex+age+severity+comorbidity,data=dev_data,family=binomial(link="logit"))
summary(reg)
#> 
#> Call:
#> glm(formula = y ~ sex + age + severity + comorbidity, family = binomial(link = "logit"), 
#>     data = dev_data)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.3801  -0.9824   0.4624   0.8912   2.2420  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.9094     0.1752  -5.192 2.08e-07 ***
#> sex           1.1167     0.2545   4.388 1.14e-05 ***
#> age           0.5434     0.1136   4.785 1.71e-06 ***
#> severity      0.4414     0.0599   7.368 1.73e-13 ***
#> comorbidity   0.8953     0.2079   4.306 1.66e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 686.86  on 499  degrees of freedom
#> Residual deviance: 564.21  on 495  degrees of freedom
#> AIC: 574.21
#> 
#> Number of Fisher Scoring iterations: 4
```

Given this, our risk prediction model can be written as:

*l**o**g**i**t*(*p*) =  − 0.9094 + 1.1167 \* *s**e**x* + 0.5434 \* *a**g**e* + 0.4414 \* *s**e**v**e**r**i**t**y* + 0.8953 \* *c**o**m**o**r**b**i**d**i**t**y*.

First, let’s compare the ROC and mROC in the development data

``` r
pred<-predict.glm(reg, type='response')

library(pROC)
#> Type 'citation("pROC")' for a citation.
#> 
#> Attaching package: 'pROC'
#> The following objects are masked from 'package:stats':
#> 
#>     cov, smooth, var
dev_roc<-roc(response=dev_data[,'y'], predictor=pred)
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
plot(dev_roc)
title("ROC in the development dataset")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

And now the much-awaited mROC using these data. Note that we use the
line function to add the mROC on top

``` r
dev_mroc<-mROC(p=pred)
```

``` r
plot(dev_roc)
lines(dev_mroc, col="red")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Important note: the statistical inference on comparing mROC and ROC cannot be used for internal validation. Such a test is made for external validation.

Now lets calculate the predicted probabilities for each subject in the
validation dataset given the prediction equation.

``` r
pred<-predict.glm(reg,newdata = val_data, type="response")

summary(pred)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02774 0.37566 0.57743 0.55688 0.73410 0.98611
```

Using the package pROC, let’s draw the validation ROC curve

``` r
val_roc<-roc(response=val_data[,'y'], predictor=pred)
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
plot(val_roc)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

And now the much-awaited mROC using these data. Note that we use the
line function to add the mROC on top

``` r
val_mroc<-mROC(p=pred)
```

Notice that the mROC function only requires the vector of predicted
probabilities.

To compare the ROC and mROC plots, we juxtapose them next to each other:

``` r
plot(val_roc)
lines(val_mroc, col="red")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Here, it is obvious that the mROC and ROC curve are not compatible,
indicating that the model is not calibrated.

``` r
res<-mROC_inference(val_data[,'y'],pred)

res
#> Mean calibration statistic (A):0.08187776(Obs<Pred) (p:0.00019)
#> mROC/ROC equality statsitic (B):0.06583836 (p:0.00371)
#> Unified statistic:28.50314 (df:4.012282,p:1.000175e-05)
```
