
hdma <- function (X, Y, M, COV.XM = NULL, COV.MY = COV.XM, family = c("gaussian","binomial"), method = c("lasso", "ridge"), topN = NULL,
                  parallel = FALSE, ncore = 1, verbose = FALSE, ...){
  pkgs <- list("hdi","MASS","doParallel", "foreach","iterators")
  checking<-unlist(lapply(pkgs, require, character.only = T))
  if(any(checking==F))
    stop("Please install the necessary packages first!")
  family <- match.arg(family)
  method <- match.arg(method)
  if (parallel & (ncore == 1)) ncore <- parallel::detectCores()
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2*n/log(n))
  } else {
    d <- topN      # the number of top mediators that associated with independent variable (X)
  }
  d <- min(p, d)   # if d > p select all mediators
  #############################################################################################################################
  ################################           Step-1 Sure Independent Screening (SIS)          #################################
  #############################################################################################################################
  message("Step 1: Sure Independent Screening ...", " (", Sys.time(), ")")
  if(family == "binomial")
  {
    if(verbose) message("Screening M using the association between X and M: ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(Y=NA, M, X, COV=COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", parallel = parallel,
                                   ncore = ncore, verbose, tag = "Sure Independent Screening")


    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", parallel = parallel,
                           ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)
  #################################################################################################################################
  ####################################          Step-2  High-dimensional Inference (HDI)         ##################################
  #################################################################################################################################
  message("Step 2: High-dimensional inference (", method, ") ...", "     (", Sys.time(), ")")
  ## Based on the SIS results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM, Y, family = family) else fit <- ridge.proj(XM, Y, family = family)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM_COV, Y, family = family) else fit <- ridge.proj(XM_COV, Y, family = family)
  }
  P_hdi<-fit$pval[1:length(ID)]
  index<-which(P_hdi<=0.05)
  if(verbose)  message("Non-zero ",method, " beta estimate(s) of mediator(s) found: ", paste0(names(index), collapse = ","))
  if(length(index)==0)
    stop("No mediatior is identified !")
  ID_test <- ID[index]
  if(family == "binomial")
  {
    ## This has been done in step 1 (when Y is binary)
    alpha_est <- alpha[,ID_test, drop = FALSE]
  } else {
    if(verbose) message(" Estimating alpha (effect of X on M): ", appendLF = FALSE)
    alpha_est <- himasis(NA, M[,ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X",
                         parallel = FALSE, ncore = ncore, verbose, tag = "site-by-site ordinary least squares estimation")
  }
  beta_P<-P_hdi[index]
  beta_hat<-fit$bhat[index]              # the estimator for beta
  alpha_hat<-as.numeric(alpha_est[1, ])
  ab_est<-beta_hat*alpha_hat
  alpha_P<-alpha_est[2,]
  PA <- rbind(beta_P,alpha_P)
  P_value <- apply(PA,2,max)
  ###################################################################################################################################
  ###############################          STEP 3   Computing the propotion of mediation effect          ############################
  ###################################################################################################################################
  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
  results <- data.frame(alpha = alpha_hat, beta = beta_hat, gamma = gamma_est, `alpha*beta` = ab_est, `%total effect`
                        =ab_est/gamma_est*100, `P.value` = P_value, check.names = FALSE)
  rownames(results) <- ID_test
  message("Done!", " (", Sys.time(), ")")
  doParallel::stopImplicitCluster()
  return(results)
}
