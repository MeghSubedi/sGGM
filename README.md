
<!-- README.md is generated from README.Rmd. Please edit that file -->

   

This package provides tools to construct sparse Phenotype-Phenotype
Network (PPN) from GWAS summary statistics. It enables researchers to
infer conditional dependencies among phenotypes by estimating sparse
precision matrices using regularized methods. This package is
particularly useful for identifying phenotype clusters,visualization of
the network structure and multiple phenotype association tests.

## Installation

You can install the development version of sGGM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("MeghSubedi/sGGM")
```

## Example

One built-in example dataset is included with the package to facilitate
exploration of the sGGM framework. This dataset contains 10,000 rows and
40 columns, representing GWAS summary statistics for 10,000 SNPs across
40 phenotypes.

``` r
# load the library
library(sGGM)

# load the example data 
data("example_data")
dim(example_data)
#> [1] 10000    40
```

``` r
example_data[1:10,1:8]
#>          V1       V2       V3       V4       V5       V6       V7       V8
#> 1  2.352686 3.102459 4.659811 3.952957 2.196122 2.544713 2.753868 3.241895
#> 2  3.972366 2.536698 4.109112 3.540926 2.002657 2.981733 4.405845 3.070740
#> 3  7.661438 4.969129 8.355886 7.091015 5.219148 5.656682 5.094139 4.469091
#> 4  5.495580 5.841046 7.393212 5.947068 3.558339 4.095446 4.319241 4.677313
#> 5  5.898467 3.805505 5.470657 5.970598 5.775413 4.521835 5.775413 5.596331
#> 6  4.625398 3.920994 3.854789 3.148123 2.240568 2.613996 2.240568 2.053854
#> 7  4.842526 3.104963 4.693197 4.002270 5.264937 3.982452 2.767467 3.239961
#> 8  4.543508 3.968040 4.901397 4.651954 2.843762 4.364845 4.629381 3.306700
#> 9  5.870082 5.366826 5.984522 5.458167 4.882080 5.574574 6.890312 5.505324
#> 10 4.998970 5.309177 5.785537 6.669334 4.625818 4.047591 4.520686 4.257855
```

Lets explore the cluster structure obtained from the correlation matrix
and partial correlation matrix. Get_Clusters() function will return the
list containing cluster structure from correlation matrix as well as
from partial correlation matrix.

``` r
# Detect clusters from correlation and partial correlation matrix 
Modules<-Get_Clusters(example_data)

# The cluster structure from the correlation matrix 
Modules$B_cor
#>       [,1] [,2] [,3] [,4]
#>  [1,]    1    0    0    0
#>  [2,]    1    0    0    0
#>  [3,]    1    0    0    0
#>  [4,]    1    0    0    0
#>  [5,]    1    0    0    0
#>  [6,]    1    0    0    0
#>  [7,]    1    0    0    0
#>  [8,]    1    0    0    0
#>  [9,]    1    0    0    0
#> [10,]    1    0    0    0
#> [11,]    0    1    0    0
#> [12,]    0    1    0    0
#> [13,]    0    1    0    0
#> [14,]    0    1    0    0
#> [15,]    0    1    0    0
#> [16,]    0    1    0    0
#> [17,]    0    1    0    0
#> [18,]    0    1    0    0
#> [19,]    0    1    0    0
#> [20,]    0    1    0    0
#> [21,]    0    0    1    0
#> [22,]    0    0    1    0
#> [23,]    0    0    1    0
#> [24,]    0    0    1    0
#> [25,]    0    0    1    0
#> [26,]    0    0    1    0
#> [27,]    0    0    1    0
#> [28,]    0    0    1    0
#> [29,]    0    0    1    0
#> [30,]    0    0    1    0
#> [31,]    0    0    0    1
#> [32,]    0    0    0    1
#> [33,]    0    0    0    1
#> [34,]    0    0    0    1
#> [35,]    0    0    0    1
#> [36,]    0    0    0    1
#> [37,]    0    0    0    1
#> [38,]    0    0    0    1
#> [39,]    0    0    0    1
#> [40,]    0    0    0    1
```

``` r
# The cluster structure from partial correlation matrix 
Modules$B_Pcor
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [5,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [6,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [7,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [8,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [9,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#> [10,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#> [11,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [12,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [13,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [14,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [15,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [16,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [17,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [18,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [20,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [21,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#> [22,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#> [23,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [24,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [25,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [26,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [27,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [28,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [29,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [30,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [31,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [32,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [33,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [34,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [35,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [36,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [37,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [38,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [39,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> [40,]    0    0    0    0    0    0    0    0    0     0     0     0     1
```

Now, we can use these network modules to test the association of SNPs
with multiple phenotype.

``` r
# Association test without considering Network modules 
Pval<-Test_Without_Network(example_data,"Wald")
# This will return vector of length 10000 equivalent to the number of SNPs 
length(Pval)
#> [1] 10000
```

``` r
head(Pval)
#>            1            2            3            4            5            6 
#> 1.560975e-03 6.589561e-02 1.049159e-09 9.416528e-07 4.270389e-06 1.574251e-02
```

``` r
# Association test considering network modules from correlation matrix 
Pval<-Test_With_Network(example_data,Modules$B_cor,method="Wald")
head(Pval)
#>       Module 1   Module 2   Module 3   Module 4
#> 1 4.699393e-03 0.01715324 0.61916979 0.03690649
#> 2 1.042088e-03 0.54056479 0.65517078 0.73515994
#> 3 2.864664e-14 0.55857769 0.22529293 0.45920470
#> 4 3.083545e-10 0.68112187 0.02223023 0.99595896
#> 5 5.704259e-09 0.45042021 0.02228034 0.55732716
#> 6 4.430775e-03 0.59114815 0.20315658 0.12799165
```

``` r
dim(Pval)
#> [1] 10000     4
```

The result is a data frame containing p-values from association tests
conducted on phenotypes within each module. These p-values can then be
combined to obtain a final p-value for the global association test.

``` r
# Association test considering network modules from partial correlation matrix 
Pval<-Test_With_Network(example_data,Modules$B_Pcor,"Wald")
head(Pval)
#>       Module 1     Module 2     Module 3   Module 4  Module 5  Module 6
#> 1 6.632954e-03 2.033332e-04 1.167537e-02 0.03401537 0.1343972 0.5150347
#> 2 1.685403e-03 5.032044e-04 5.402343e-05 0.88188860 0.1697812 0.4862409
#> 3 1.409980e-12 1.805985e-14 1.017301e-08 0.94925424 0.8788995 0.1198473
#> 4 2.428318e-08 4.775480e-13 1.384670e-05 0.94899000 0.1432924 0.4088750
#> 5 2.078626e-10 3.752407e-09 2.143730e-08 0.94406133 0.3476275 0.4584150
#> 6 1.330541e-04 2.643354e-03 2.242036e-02 0.26921498 0.8680351 0.1943001
#>    Module 7    Module 8  Module 9  Module 10 Module 11 Module 12  Module 13
#> 1 0.2906282 0.787623990 0.8222264 0.01278614 0.3556652 0.5945651 0.02855782
#> 2 0.6975374 0.272777032 0.5860767 0.99513323 0.4336538 0.2773428 0.97536022
#> 3 0.2636178 0.229189514 0.1331612 0.11960471 0.3219397 0.5624201 0.55887030
#> 4 0.3576047 0.006085071 0.2267573 0.98582571 0.8020317 0.6894035 0.97843492
#> 5 0.3001211 0.016609932 0.1194859 0.62808166 0.6947809 0.5675131 0.88107900
#> 6 0.2830739 0.335859987 0.1115829 0.13296817 0.3640034 0.8934061 0.33392377
```
