# qtQDA: quantile transformed quadratic discriminant analysis for high-dimensional RNA-seq data
qtQDA contains functions for classification, quantile transformation and quadratic discrimant analysis and several generalized linear model analysis of RNA-Seq data. 

## Installation

``` r
# install.packages("devtools")
devtools::install_github("goknurginer/qtQDA")
library(qtQDA)
```

For complete list of functions and instructions: 

```r
library(help = "qtQDA") 
```
Note that devtools does not build vignettes by default. To view the vignette: 
```r
devtools::install_github("goknurginer/qtQDA", build_vignettes = TRUE)
library(qtQDA)
vignette("qtQDA")
```

# Citation:
The article explaining the method can be found and cited [here](https://peerj.com/articles/8260/).
