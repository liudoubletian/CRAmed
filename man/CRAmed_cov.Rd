% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRAmed_cov.R
\name{CRAmed_cov}
\alias{CRAmed_cov}
\title{CRAmed: A conditional randomization test for sparse and high-dimensional mediation analysis in microbiome data}
\usage{
CRAmed_cov(
  M_mat,
  Y,
  Exposure,
  FDR = 0.05,
  n.times = 100,
  prefilter = TRUE,
  cov,
  n.perm = 100,
  CI = FALSE,
  modely = "gaussian",
  modelm = "ZINB",
  method = "BH"
)
}
\arguments{
\item{M_mat}{a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.}

\item{Y}{a continuous outcome.}

\item{Exposure}{a binary exposure.}

\item{FDR}{the threshold value for p-values, with a default of 0.05.}

\item{n.times}{the number of times for sampling mediators from the ZINB model, with a default of 100.}

\item{prefilter}{If TRUE, the lasso method will be used to select variables; default is TRUE.}

\item{cov}{the matrix containing covariates.}

\item{n.perm}{the number of times for sampling dataset for calculating confidence interval, with a default of 100.}

\item{CI}{If TRUE, the confidence intervals for NDE, NIE, NIEA, and NIEP will be calculated; default is FALSE.}

\item{modely}{The model for fitting the outcome variable Y, including 'gaussian' or 'binomial'; default is 'gaussian'.}

\item{modelm}{The model for fitting the mediator M, including 'ZINB', 'NB', and 'ZIP'; default is 'ZINB'.}

\item{method}{The method for correcting for multiple testing, with a default of 'BH' (Benjamini-Hochberg).}
}
\value{
Return the identified mediators, corresponding p-values and confidence interval estimates of NDE, NIE, NIEA, NIEP, and AIC, BIC, residuals for the corresponding models.
}
\description{
CRAmed: A conditional randomization test for sparse and high-dimensional mediation analysis in microbiome data
}
\examples{
otu_n <- 50;num <- 50
set.seed(123)
sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)
cramedcov.res <- CRAmed_cov(M_mat=sim_zinb.mat$M_mat,Y=sim_zinb.mat$Y, Exposure=sim_zinb.mat$trt, cov=matrix(rnorm(num)), n.perm=10, CI=TRUE)
}
