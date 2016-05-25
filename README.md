
# condSURV: Estimation of the Conditional Survival Function for Ordered Multivariate Failure Time Data

```condSURV``` is an R package to implement some newly developed methods 
for the estimation of the conditional survival function. The package implements
three nonparametric and semiparametric estimators for these quantities. The package also implements feasible estimation methods for these quantities conditionally on current or past covariate measures. 

Other related estimators are also implemented in the package. One of these estimators is the Kaplan-Meier estimator typically assumed to estimate the survival function. A modification of the Kaplan-Meier
estimator, based on a preliminary estimation (presmoothing) of the
censoring probability for the survival time, given the available information is also implemented.


## Installation
```condSURV``` is available through both CRAN and GitHub.

Get the released version from CRAN:
```
install.packages("condSURV")
```

Or the development version from GitHub:
```
# install.packages("devtools")
devtools::install_github("sestelo/condsurv")
```
