% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Naive-CRAmed.R
\name{Creat_condition_gaussian}
\alias{Creat_condition_gaussian}
\title{Naive.CRAmed: A conditional randomization test for normal microbiome data}
\usage{
Creat_condition_gaussian(X, indx, Sigma = NULL, lambda = NULL)
}
\arguments{
\item{M_mat}{a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.}

\item{Y}{a continuous outcome.}

\item{Exposure}{a binary exposure.}

\item{FDR}{the threshold value for p-values.}

\item{method}{the correction method for multiple testing.}
}
\value{
Return the identified mediators and corresponding p-values of NIE.
}
\description{
Naive.CRAmed: A conditional randomization test for normal microbiome data
}
\examples{
otu_n <- 50;num <- 50
set.seed(12)
sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)
naive.cramed.res <- Naive.CRAmed(M_mat=log(sim_zinb.mat$M_mat+1), Y=sim_zinb.mat$Y, Exposure=sim_zinb.mat$trt)
}
