# CRAmed: A conditional randomization test for high-dimentional mediation analysis in sparse microbiome data

This is an R package from "CRAmed: A conditional randomization test for high-dimentional mediation analysis in sparse microbiome data". 

## Installation 

Before you install the CRAmed, some other libraries are required to install, see the ```requirement.R```:
  ```r
library(MASS)
library(plyr)
library(glmnet)
library(pscl)
```
Then you can install CRAmed from github with:
  ```r
install.packages("devtools")  
devtools::install_github("liudoubletian/CRAmed") 
library(CRAmed)  
```
## Vignette

The details of the manual and the code used to generate each figure in the manuscript are provided in the fold ```../vignettes```. 

Here, we show a brief example.

## Example

```r
library(CRAmed); packageVersion("CRAmed")

#Simulate the ZINB data
otu_n <- 50;num <- 50
set.seed(1)
sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)

#Detect the mediators by CRAmed
cramed.res <- CRAmed(M_mat=sim_zinb.mat$M_mat, Y=sim_zinb.mat$Y, Exposure=sim_zinb.mat$trt, n.perm=10, CI=TRUE)
cramed.res
```


