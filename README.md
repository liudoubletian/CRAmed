# CRAmed: A conditional randomization test for sparse and high-dimentional mediation analysis in microbiome data 

This is an R package from "CRAmed: A conditional randomization test for sparse and high-dimentional mediation analysis in microbiome data ". 

## Installation 

Before you install the CRAmed, some other libraries are required to install, see the ```requirement.R```:
```r
library(MASS)
library(plyr)
library(glmnet)
library(pscl)
```
Then you can also install CRAmed from github with:
```r
install.packages("devtools")  
devtools::install_github("liudoubletian/CRAmed") 
library(CRAmed)  
```
## Vignette
You can find the vignette at 


Here, we show a brief example.

Simulate an example


```r

#CRAmed: A conditional randomization test for sparse and high-dimentional mediation analysis in microbiome data 
library(CRAmed); packageVersion("CRAmed")

#Simulate the ZINB data
otu_n <- 50;num <- 50
set.seed(1)
sim_zinb.mat <- sim_zinb(otu_n,num)

#Detect the mediators by CRAmed
cramed.res <- dCRT_causal(M_mat=sim_zinb.mat$M_mat,Y=sim_zinb.mat$Y_mat, Exposure=sim_zinb.mat$t)
cramed.res
```


