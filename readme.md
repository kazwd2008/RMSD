# RMSD

<!-- badges: start -->
<!-- badges: end -->

The package RMSD contains a function RMSD. It aims to provide a tool for multivariate outlier detection for elliptical data distribution based on the Modified Stahel-Donoho (MSD) estimators based on Béguin and Hulliger(2003).

The method requires a large number of projection as the data dimension increases, so please try the paralleled version, RMSDp::RNSDp() if the function RMSD() is terminated due to insufficient memory.

The function RMSD::RMSD() is faster than the function RMSDP::RMSDp for lower dimensional datasets.

## Arguments

The parameter 'nb' controls the number of basis for projection.  Default setting is $\exp(2.1328+0.8023p)$ where $p$ is the number of variables. Specifying a smaller number will reduce memory usage while it may degrade the performance of outlier detection.

The parameter 'sd' is for a random seed. Plural execution without setting the random seed with a same dataset may yield different results, since the function RMSD() uses random number.

The default threshold set to decide outlier is 99.9 percentile point of F-statistics. It can be changed using the parameter 'pt'.

This function is an improved version of 'msd' at https://github.com/kazwd2008/MSD by adding the last step to decide outliers by calculating the squared Mahalanobis distance of each observation and F-statistics to decide outliers.

## Value

- u	  : estimated robust mean vector 
- V	  : estimated robust covariance matrix 
- wt  : final weights 
- mah : squared Mahalanobis distnace of each observation
- FF  : F test statistics
- cf  : threshold to detect outliers (percentile point) 
- ot  : outlier flag (1:normal observation, 2:outlier)

## Installation

### from CRAN

The package RMSDP is on CRAN.

``` r
install.packages("RMSDps")
library(RMSDp)
```

### from github

You can also install the package from [GitHub](https://github.com/kazwd2008) with:

``` r
install.packages("devtools")
devtools::install_github("kazwd2008/RMSDp")
```

## Example: wine data [13 variables]

This is an example with the bushfire data set used by Campbell(1984) to locate bushfire scars.

The dataset contains 38 observations with 5 variables.

## Usage

```
# install.packages("RMSD")
require(RMSD)

data(bushfire, package="robustbase")	
o.msd <- RMSD(as.matrix(bushfire)) 

# scatterplot matrix
pairs(bushfire, pch=19, col=o.msd$ot)

# Q-Q plot
n <- nrow(bushfire)
d <- ncol(bushfire)
qqplot(qchisq(ppoints(n), df=d), o.msd$mah, pch=19, col=sort(o.msd$ot), main = "Q-Q plot")
abline(0, 1, col = 'green')
```

## Reference

Béguin, C. and B. Hulliger (2003) Robust Multivariate Outlier Detection and Imputation with Incomplete Survey Data, EUREDIT Deliverable D4/5.2.1/2 Part C. 

Campbell, N. A.（1989）, Bushfire mapping using noaa avhrr data, Technical report, CSIRO
 