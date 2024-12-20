
#' CRAmed: A conditional randomization test for sparse and high-dimensional mediation analysis in microbiome data
#'
#' @param M_mat a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param Y a continuous outcome.
#' @param Exposure a binary exposure.
#' @param cov the matrix containing covariates.
#' @param FDR the threshold value for p-values, with a default of 0.05.
#' @param n.times the number of times for sampling mediators from the ZINB model, with a default of 100.
#' @param n.perm  the number of times for sampling dataset for calculating confidence interval, with a default of 100.
#' @param CI If TRUE, the confidence intervals for NDE, NIE, NIEA, and NIEP will be calculated; default is FALSE.
#' @param prefilter If TRUE, the lasso method will be used to select variables; default is TRUE.
#' @param modely The model for fitting the outcome variable Y, including 'gaussian' or 'binomial'; default is 'gaussian'.
#' @param modelm The model for fitting the mediator M, including 'ZINB', 'NB', and 'ZIP'; default is 'ZINB'.
#' @param method The method for correcting for multiple testing, with a default of 'BH' (Benjamini-Hochberg).

#' @return Return the identified mediators, corresponding p-values and confidence interval estimates of NDE, NIE, NIEA, NIEP, and AIC, BIC, residuals for the corresponding models.
#' @examples
#' otu_n <- 50;num <- 50
#' set.seed(123)
#' sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)
#' cramedcov.res <- CRAmed_cov(M_mat=sim_zinb.mat$M_mat,Y=sim_zinb.mat$Y, Exposure=sim_zinb.mat$trt, cov=matrix(rnorm(num)), n.perm=10, CI=TRUE)
#' @export
#' @importFrom glmnet  cv.glmnet
#' @import MASS
#' @import plyr
#' @import pscl


CRAmed_cov <- function(M_mat, Y, Exposure, FDR = 0.05, n.times=100, prefilter=TRUE, cov,
                   n.perm=100, CI=FALSE, modely = 'gaussian', modelm = 'ZINB', method = 'BH') {
  library(glmnet)
  library(MASS)
  Creat_condition_zinb_cov <- function(M_mat, indx, Exposure, n.times, cov){
    library(pscl)
    p <- ncol(M_mat)
    n <- nrow(M_mat)
    coe_mat <- c()
    Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(M_mat[,indx]),cov=cov)
    zinb.fit <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, data = Y_data,dist = c("negbin"),link = c("logit")), silent = TRUE)
    if(inherits(zinb.fit,"try-error")){
      zinb.fit <-  try(glm.nb(Mediator ~ Exposure+cov, data = Y_data), silent = TRUE)
      alpha <- summary(zinb.fit)$coefficients[,1]
      phi <- zinb.fit$theta

      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      x <- cbind(1,Exposure,cov)
      lamda <- exp(x%*%matrix(alpha))
      samp_M <- matrix(NA, nrow(M_mat), n.times)
      for(j in 1:n.times) {
        set.seed(j)
        for (i in 1:n) {
          samp_M[i,j] = rnbinom(1, mu=lamda[i], size=phi)
        }
      }
      res_M_samp <- lapply(c(1:n.times), function(j){
        Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(samp_M[,j]),cov=cov)
        zinb.fit.samp <-  try(glm.nb(Mediator ~ Exposure+cov, data = Y_data), silent = TRUE)
        if(inherits(zinb.fit.samp,"try-error")){j <- j+1 }
        else{data.frame(t(zinb.fit.samp$residuals))}
      })
    }
    else{
      alpha <- summary(zinb.fit)$coefficients$count[-(ncol(Y_data)+1),1]
      gamma <- summary(zinb.fit)$coefficients$zero[-(ncol(Y_data)+1),1]
      phi <- zinb.fit$theta

      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      x <- cbind(1,Exposure,cov)
      logit.p  <- x %*% matrix(gamma)
      p  <- 1 / (1 + exp(-logit.p))
      lamda <- exp(x%*%matrix(alpha))
      samp_M <- matrix(NA, nrow(M_mat), n.times)
      for(j in 1:n.times) {
        set.seed(j)
        Z1 <- rbinom(n, 1, p)
        for (i in 1:n) {
          if (Z1[i] == 1) {
            samp_M[i,j] = 0
          }
          if (Z1[i] == 0) {
            samp_M[i,j] = rnbinom(1, mu=lamda[i], size=phi)
          }
        }
      }
      res_M_samp <- lapply(c(1:n.times), function(j){
        if(sum(samp_M[,j])==0){j=j+1}
        else{
          Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(samp_M[,j]),cov=cov)
          zinb.fit.samp <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, data = Y_data,dist = c("negbin"),link = c("logit")), silent = TRUE)
          if(inherits(zinb.fit.samp,"try-error")){
            zinb.fit.samp <-  try(glm.nb(Mediator ~ Exposure+cov, data = Y_data), silent = TRUE)
            if(inherits(zinb.fit.samp,"try-error")){
              zinb.fit.samp <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, data = Y_data,dist = c("poisson"),link = c("logit")), silent = TRUE)
            }
          }
          data.frame(t(zinb.fit.samp$residuals))}
      }
      )

    }

    return(list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp))
  }
  Creat_condition_zip_cov <- function(M_mat, indx, Exposure, n.times, cov){
    library(pscl)
    p <- ncol(M_mat)
    n <- nrow(M_mat)
    coe_mat <- c()
    Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(M_mat[,indx]),cov=cov)
    zip.fit <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, data = Y_data,dist = c("poisson"),link = c("logit")), silent = TRUE)

      alpha <- summary(zip.fit)$coefficients$count[-(ncol(Y_data)+1),1]
      gamma <- summary(zip.fit)$coefficients$zero[-(ncol(Y_data)+1),1]

      M_bar <- zip.fit$fitted.values
      M_res <- zip.fit$residuals
      x <- cbind(1,Exposure,cov)
      logit.p  <- x %*% matrix(gamma)
      p  <- 1 / (1 + exp(-logit.p))
      lamda <- exp(x%*%matrix(alpha))
      samp_M <- matrix(NA, nrow(M_mat), n.times)
      for(j in 1:n.times) {
        set.seed(j)
        Z1 <- rbinom(n, 1, p)
        for (i in 1:n) {
          if (Z1[i] == 1) {
            samp_M[i,j] = 0
          }
          if (Z1[i] == 0) {
            samp_M[i,j] = rpois(1, lambda=lamda[i])
          }
        }
      }
      res_M_samp <- lapply(c(1:n.times), function(j){
        if(sum(samp_M[,j])==0){j=j+1}
        else{
          Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(samp_M[,j]),cov=cov)
          zip.fit.samp <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, data = Y_data,dist = c("poisson"),link = c("logit")), silent = TRUE)
          data.frame(t(zip.fit.samp$residuals))}
      }
      )


    return(list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp))
  }

  Creat_condition_nb_cov <- function(M_mat, indx, Exposure, n.times, cov){
    library(pscl)
    p <- ncol(M_mat)
    n <- nrow(M_mat)
    coe_mat <- c()
    Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(M_mat[,indx]),cov=cov)
      nb.fit <-  try(glm.nb(Mediator ~ Exposure+cov, data = Y_data), silent = TRUE)
      alpha <- summary(nb.fit)$coefficients[,1]
      phi <- nb.fit$theta

      M_bar <- nb.fit$fitted.values
      M_res <- nb.fit$residuals
      x <- cbind(1,Exposure,cov)
      lamda <- exp(x%*%matrix(alpha))
      samp_M <- matrix(NA, nrow(M_mat), n.times)
      for(j in 1:n.times) {
        set.seed(j)
        for (i in 1:n) {
          samp_M[i,j] = rnbinom(1, mu=lamda[i], size=phi)
        }
      }
      res_M_samp <- lapply(c(1:n.times), function(j){
        Y_data <- data.frame(Exposure=Exposure, Mediator=as.vector(samp_M[,j]),cov=cov)
        nb.fit <-  try(glm.nb(Mediator ~ Exposure+cov, data = Y_data), silent = TRUE)
        if(inherits(nb.fit,"try-error")){j <- j+1 }
        else{data.frame(t(nb.fit$residuals))}
      })


    return(list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp))
  }

  p <- ncol(M_mat)
  n <- nrow(M_mat)

  ############### CDR ###############

  if (prefilter == TRUE){
    set.seed(123)
    cv_lasso <- cv.glmnet(cbind(Exposure,cov,M_mat), Y, alpha = 1, family = modely,
                          dfmax = as.integer(p / 2))
    lamb <- cv_lasso$lambda.min
    opt_model <- glmnet(cbind(Exposure,cov,M_mat), Y, alpha = 1, lambda = lamb,
                        family = modely, dfmax = as.integer(p / 2))
    residualsy <- Y - predict(opt_model, cbind(M_mat, cov, Exposure))
    beta_fit <- opt_model$beta[-(1:(ncol(cov)+1))]
    beta_2 <- opt_model$beta[1]
    selection_set <- which(beta_fit != 0)
  }else{
    y.data <- data.frame(Y=Y,Exposure,cov,M_mat)
    lm.fit <- lm(Y~.,data=y.data)
    residualsy <- lm.fit$residuals
    beta_fit <- summary(lm.fit)$coefficients[-(1:(ncol(cov)+1)),1]
    selection_set <- 1:p
  }


  ############## dCRT ##############
  pvl_lst <- c()
  est_lst <- c()
  if(length(selection_set)==0){pvl_lst=NA}
  else{
    if(p==1){
      if(modelm == 'ZINB'){Cond_M <- Creat_condition_zinb_cov(M_mat, indx=1, Exposure, n.times, cov=cov)}
      if(modelm == 'ZIP'){Cond_M <- Creat_condition_zip_cov(M_mat, indx=1, Exposure, n.times, cov=cov)}
      if(modelm == 'NB'){Cond_M <- Creat_condition_nb_cov(M_mat, indx=1, Exposure, n.times, cov=cov)}
      mean_M <- Cond_M$mean_m
      M_res_ob <- Cond_M$res_m
      data <- data.frame(Y=Y,Exposure=Exposure,cov=cov)
      if (modely == 'binomial'){
        eps_res <- Y - 1 / (1 + exp(- predict(lm(Y~Exposure+cov,family = modely,data=data))))
      }
      if (modely == 'gaussian'){
        eps_res <- Y - predict(lm(Y~Exposure+cov,family = modely,data=data))
      }
      imp_obe <- abs(mean(M_res_ob * eps_res)) / mean(M_res_ob^2)
      ############### Estimate p-values ##################
      library(plyr)
      list.s <- unlist(lapply(Cond_M$res_m_samp,function(x){length(x)}))
      list.sel <-which(list.s==1)
      if(length(list.sel)!=0){ M_res_sample <- as.matrix(do.call(rbind.fill,Cond_M$res_m_samp[-list.sel]))
      }
      else{
        M_res_sample <- as.matrix(do.call(rbind.fill,Cond_M$res_m_samp))
      }
      var_lst_sample <- unlist(lapply(c(1:nrow(M_res_sample)), function(j){
        mean((unlist(M_res_sample[j,]))^2) }))
      t_lst <- abs(M_res_sample %*% eps_res / n) / var_lst_sample
      pvl_lst[1] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
      beta_fit <- NULL

      print('############# d0CRT pvalue ###############')
      print(paste(1, pvl_lst[1], sep = ': '))
    }
    ##########
    if(p>1){

      Cond_E <- list()
      for(j in 1:n.times){
        set.seed(j)
        indx.e <- sample(1:n)
        Cond_E$res_e_samp[[j]] <- as.data.frame(t(Exposure[indx.e]))
      }
      ############### Distill Y and Importance ##################
      set.seed(123)
      cv_lasso_null <- cv.glmnet(cbind(M_mat,cov), Y, alpha = 1, family = modely,
                                 dfmax = as.integer(p / 2))
      lamb_null <- cv_lasso_null$lambda.min
      model_res_null <- glmnet(cbind(M_mat,cov), Y, alpha = 1, lambda = lamb_null,
                               family = modely, dfmax = as.integer(p / 2))
      if (modely == 'binomial'){
        eps_res <- Y - 1 / (1 + exp(- predict(model_res_null,cbind(M_mat,cov))))
      }
      if (modely == 'gaussian'){
        eps_res <- Y - predict(model_res_null, cbind(M_mat,cov))
      }

      imp_exp <- abs(mean(Exposure * eps_res)) / mean(Exposure^2)

      ############### Estimate p-values ##################
      library(plyr)
      list.s <- unlist(lapply(Cond_E$res_e_samp,function(x){length(x)}))
      list.sel <-which(list.s==1)
      if(length(list.sel)!=0){ 
        E_res_sample <- as.matrix(do.call(rbind.fill,Cond_E$res_e_samp[-list.sel]))
      }
      else{
        E_res_sample <- as.matrix(do.call(rbind.fill,Cond_E$res_e_samp))
      }
      var_lst_sample <- unlist(lapply(c(1:nrow(E_res_sample)), function(j){
        mean((unlist(E_res_sample[j,]))^2) }))
      t_lst <- abs(E_res_sample %*% eps_res / n) / var_lst_sample
      nde.p <- mean(c(1, ifelse(t_lst >= imp_exp, 1, 0)))

      for (j in 1:length(selection_set)){

        indx <- selection_set[j]

        ############### Distill X ##################
        if(modelm == 'ZINB'){Cond_M <- Creat_condition_zinb_cov(M_mat, indx, Exposure, n.times, cov=cov)}
        if(modelm == 'ZIP'){Cond_M <- Creat_condition_zip_cov(M_mat, indx, Exposure, n.times, cov=cov)}
        if(modelm == 'NB'){Cond_M <- Creat_condition_nb_cov(M_mat, indx, Exposure, n.times, cov=cov)}
        mean_M <- Cond_M$mean_m
        M_res_ob <- Cond_M$res_m

        ############### Distill Y and Importance ##################
        set.seed(123)
        cv_lasso_null <- cv.glmnet(cbind(M_mat[,-indx],Exposure,cov), Y, alpha = 1, family = modely,
                                   dfmax = as.integer(p / 2))
        lamb_null <- cv_lasso_null$lambda.min
        model_res_null <- glmnet(cbind(M_mat[,-indx],Exposure,cov), Y, alpha = 1, lambda = lamb_null,
                                 family = modely, dfmax = as.integer(p / 2))
        if (modely == 'binomial'){
          eps_res <- Y - 1 / (1 + exp(- predict(model_res_null,cbind(M_mat[,-indx],Exposure,cov))))
        }
        if (modely == 'gaussian'){
          eps_res <- Y - predict(model_res_null, cbind(M_mat[,-indx],Exposure,cov))
        }

        imp_obe <- abs(mean(M_res_ob * eps_res)) / mean(M_res_ob^2)
        #est_lst[indx] <-

        ############### Estimate p-values ##################
        library(plyr)
        list.s <- unlist(lapply(Cond_M$res_m_samp,function(x){length(x)}))
        list.sel <-which(list.s==1)
        if(length(list.sel)!=0){ M_res_sample <- as.matrix(do.call(rbind.fill,Cond_M$res_m_samp[-list.sel]))
        }
        else{
          M_res_sample <- as.matrix(do.call(rbind.fill,Cond_M$res_m_samp))
        }
        var_lst_sample <- unlist(lapply(c(1:nrow(M_res_sample)), function(j){
          mean((unlist(M_res_sample[j,]))^2) }))
        t_lst <- abs(M_res_sample %*% eps_res / n) / var_lst_sample
        pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))

        print('############# d0CRT pvalue ###############')
        print(paste(indx, pvl_lst[indx], sep = ': '))

      }
    }

  tij_mat <- log(rowSums(M_mat))
  index.p <- pvl_lst[selection_set]
  final.p <- c()
  aic <- c()
  bic <- c()
  niea <- niep <- nie <- c()
  alpha.p <- gamma.p <- c()
  residualsm <- list()
  alpha_mat <-  gamma_mat <- matrix(NA,length(selection_set),2+ncol(cov))
  if(modelm=="ZINB"){
  for(j in 1:length(selection_set)){
    library(MASS)
    Y_data <- data.frame(Exposure=Exposure,Mediator=as.vector(M_mat[,selection_set[j]]),cov=cov)

    if(sum(Y_data$Mediator==0)==0){
      nb.fit <-  glm.nb(Mediator ~ Exposure+cov, data = Y_data)
      residualsm[[j]] <- nb.fit$residuals
      wald.p <- summary(nb.fit)$coefficients[2,4]
      final.p[j] <- wald.p
    }
    else{
      zinb.fit <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, offset=tij_mat,data = Y_data,dist = c("negbin"),link = c("logit")), silent = TRUE)
      if(inherits(zinb.fit, "try-error")){final.p[j] <- NA}
      else{
        aic[j] <-  AIC(zinb.fit); bic[j] <-  BIC(zinb.fit)
        residualsm[[j]] <- zinb.fit$residuals
        alpha_mat[j,]<- summary(zinb.fit)$coefficients$count[1:(ncol(cov)+2),1]
        gamma_mat[j,] <- summary(zinb.fit)$coefficients$zero[1:(ncol(cov)+2),1]
        alpha.p[j] <- summary(zinb.fit)$coefficients$count[2,4]
        gamma.p[j] <- summary(zinb.fit)$coefficients$zero[2,4]
        cov.mat<- matrix(c(vcov(zinb.fit)[2,2],vcov(zinb.fit)[2,4+ncol(cov)],vcov(zinb.fit)[4+ncol(cov),2],vcov(zinb.fit)[4+ncol(cov),4+ncol(cov)]),2,2)
        wald.t <- try(t(matrix(c(alpha_mat[j,2],gamma_mat[j,2])))%*%ginv(cov.mat)%*%matrix(c(alpha_mat[j,2],gamma_mat[j,2])), silent = TRUE)
        if(inherits(wald.t, "try-error")){final.p[j] <- NA;}
        else{
          wald.p <- 1-pchisq(as.numeric(wald.t),df=2)
          final.p[j] <- wald.p
        }

      }

    }

    beta.val <- beta_fit[selection_set[j]];
    alpha.val <- alpha_mat[j,];
    gamma.val <- gamma_mat[j,];
    niea[j] <- mean(beta.val*(1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val[-2]))))*
                      (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val)+tij_mat)-
                         exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val[-2])+tij_mat)))
    niep[j] <- mean(beta.val*(exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val)+tij_mat))*
                      ((1/(1+exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val))))-
                         (1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val[-2]))))))

    nie[j] <-niea[j]+niep[j]

  }
  bind.p <- rbind(p.adjust(index.p,method),p.adjust(final.p,method))
  joint.p <- apply(bind.p, 2, max)
  index.mi <- which(joint.p <= FDR)
  mediation_set <- selection_set[index.mi]

  nie <- nie[index.mi]
  niea <- niea[index.mi]
  niep <- niep[index.mi]

  nie.p <- joint.p[index.mi]

  bind.p <- rbind(p.adjust(index.p,method),p.adjust(alpha.p,method))
  joint.p <- apply(bind.p, 2, max)
  niea.p <- joint.p[index.mi]

  bind.p <- rbind(p.adjust(index.p,method),p.adjust(gamma.p,method))
  joint.p <- apply(bind.p, 2, max)
  niep.p <- joint.p[index.mi]
  if(CI==TRUE){
    if(length(mediation_set)==0){
      nde.perm <- c()
      Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)
      for(jk in 1:n.perm){
        set.seed(jk+2000)
        Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
        tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
        M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
        t.perm <- Y_data.perm$Exposure
        Y_mat.perm <- Y_data.perm$Outcom
        cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]
        set.seed(123)#set.seed(123)
        cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                              lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
        lamb <- cv_lasso$lambda.min
        opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                            family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
        bet.vec <- opt_model$beta[mediation_set]
        nde.perm[jk] <- opt_model$beta[(1+p)]

        print(jk)
      }
      ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

      nie = niea = niep = nie.p = niea.p = niep.p = ci.nie = ci.niea = ci.niep = NA}else{
        #######confidence interval ######
        nde.perm <- c()
        niea.perm <- matrix(NA,n.perm,length(mediation_set))
        niep.perm <- matrix(NA,n.perm,length(mediation_set))
        nie.perm <- matrix(NA,n.perm,length(mediation_set))
        Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)

        for(jk in 1:n.perm){
          set.seed(jk+2000)
          Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
          tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
          M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
          t.perm <- Y_data.perm$Exposure
          Y_mat.perm <- Y_data.perm$Outcom
          cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]

          set.seed(123)#set.seed(123)
          cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                                lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
          lamb <- cv_lasso$lambda.min
          opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                              family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
          bet.vec <- opt_model$beta[mediation_set]
          nde.perm[jk] <- opt_model$beta[(1+p)]

          for(j in 1:length(mediation_set)){
            Y_data.dcrt <- as.data.frame(cbind(m.fi = as.vector(M_mat.perm[,mediation_set[j]]),Exposure=t.perm,Outcom=Y_mat.perm,Cov=cov.perm))

            zinb.fit <-  try(zeroinfl(m.fi ~ Exposure+Cov | Exposure+Cov,offset=tij_mat.perm, data = Y_data.dcrt,dist = c("negbin"),link = c("logit")), silent = TRUE)
            if(inherits(zinb.fit,"try-error")){niea.perm[jk,j]=niep.perm[jk,j]=nie.perm[jk,j]=NA}
            else{
              beta.val <- bet.vec[j]
              alpha.val<- summary(zinb.fit)$coefficients$count[1:(2+ncol(cov)),1]
              gamma.val <- summary(zinb.fit)$coefficients$zero[1:(2+ncol(cov)),1]


              niea.perm[jk,j] <- mean(beta.val*(1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(gamma.val[-2]))))*
                                        (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(alpha.val)+tij_mat.perm)-
                                           exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(alpha.val[-2])+tij_mat.perm)))
              niep.perm[jk,j] <- mean(beta.val*(exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(alpha.val)+tij_mat.perm))*
                                        ((1/(1+exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(gamma.val))))-
                                           (1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(gamma.val[-2]))))))

              nie.perm[jk,j] <-niea.perm[jk,j]+niep.perm[jk,j]

            }

          }

          print(jk)
        }


        if(length(mediation_set)==1){
          ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)
          ci.niea <- quantile(as.vector(niea.perm), probs = c(0.025, 0.975),na.rm=TRUE)
          ci.niep <- quantile(as.vector(niep.perm), probs = c(0.025, 0.975),na.rm=TRUE)
          ci.nie <- quantile(as.vector(nie.perm), probs = c(0.025, 0.975),na.rm=TRUE)
        }else{
          ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

          ci.niea <- apply(niea.perm,2,function(x){
            quantile(x, probs = c(0.025, 0.975))
          })
          ci.niep <- apply(niep.perm,2,function(x){
            quantile(x, probs = c(0.025, 0.975))
          })
          ci.nie <- apply(nie.perm,2,function(x){
            quantile(x, probs = c(0.025, 0.975))
          })
        }

      }
    return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie, NIEA = niea, NIEP = niep,
                NDE.pval = nde.p, NIE.pval = nie.p, NIEA.pval = niea.p, NIEP.pval = niep.p,
                NDE.CI = ci.nde, NIE.CI = ci.nie, NIEA.CI = ci.niea, NIEP.CI = ci.niep,
                AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
    ))
  }else{
    return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie, NIEA = niea, NIEP = niep,
                NDE.pval = nde.p, NIE.pval = nie.p, NIEA.pval = niea.p, NIEP.pval = niep.p,
                AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
    ))
  }
  }

  if(modelm=="ZIP"){
    for(j in 1:length(selection_set)){
      library(MASS)
      Y_data <- data.frame(Exposure=Exposure,Mediator=as.vector(M_mat[,selection_set[j]]),cov=cov)

      if(sum(Y_data$Mediator==0)==0){
        nb.fit <-  glm.nb(Mediator ~ Exposure+cov, data = Y_data)
        residualsm[[j]] <- nb.fit$residuals
        wald.p <- summary(nb.fit)$coefficients[2,4]
        final.p[j] <- wald.p
      }
      else{
        zip.fit <-  try(zeroinfl(Mediator ~ Exposure+cov | Exposure+cov, offset=tij_mat,data = Y_data,dist = c("poisson"),link = c("logit")), silent = TRUE)
        if(inherits(zip.fit, "try-error")){final.p[j] <- NA}
        else{
          aic[j] <-  AIC(zip.fit); bic[j] <-  BIC(zip.fit)
          residualsm[[j]] <- zip.fit$residuals
          alpha_mat[j,]<- summary(zip.fit)$coefficients$count[1:(ncol(cov)+2),1]
          gamma_mat[j,] <- summary(zip.fit)$coefficients$zero[1:(ncol(cov)+2),1]
          alpha.p[j] <- summary(zip.fit)$coefficients$count[2,4]
          gamma.p[j] <- summary(zip.fit)$coefficients$zero[2,4]
          cov.mat<- matrix(c(vcov(zip.fit)[2,2],vcov(zip.fit)[2,4+ncol(cov)],vcov(zip.fit)[4+ncol(cov),2],vcov(zip.fit)[4+ncol(cov),4+ncol(cov)]),2,2)
          wald.t <- try(t(matrix(c(alpha_mat[j,2],gamma_mat[j,2])))%*%ginv(cov.mat)%*%matrix(c(alpha_mat[j,2],gamma_mat[j,2])), silent = TRUE)
          if(inherits(wald.t, "try-error")){final.p[j] <- NA;}
          else{
            wald.p <- 1-pchisq(as.numeric(wald.t),df=2)
            final.p[j] <- wald.p
          }

        }

      }

      beta.val <- beta_fit[selection_set[j]];
      alpha.val <- alpha_mat[j,];
      gamma.val <- gamma_mat[j,];
      niea[j] <- mean(beta.val*(1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val[-2]))))*
                        (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val)+tij_mat)-
                           exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val[-2])+tij_mat)))
      niep[j] <- mean(beta.val*(exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val)+tij_mat))*
                        ((1/(1+exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val))))-
                           (1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(gamma.val[-2]))))))

      nie[j] <-niea[j]+niep[j]

    }
    bind.p <- rbind(p.adjust(index.p,method),p.adjust(final.p,method))
    joint.p <- apply(bind.p, 2, max)
    index.mi <- which(joint.p <= FDR)
    mediation_set <- selection_set[index.mi]

    nie <- nie[index.mi]
    niea <- niea[index.mi]
    niep <- niep[index.mi]

    nie.p <- joint.p[index.mi]

    bind.p <- rbind(p.adjust(index.p,method),p.adjust(alpha.p,method))
    joint.p <- apply(bind.p, 2, max)
    niea.p <- joint.p[index.mi]

    bind.p <- rbind(p.adjust(index.p,method),p.adjust(gamma.p,method))
    joint.p <- apply(bind.p, 2, max)
    niep.p <- joint.p[index.mi]
    if(CI==TRUE){
      if(length(mediation_set)==0){
        nde.perm <- c()
        Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)
        for(jk in 1:n.perm){
          set.seed(jk+2000)
          Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
          tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
          M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
          t.perm <- Y_data.perm$Exposure
          Y_mat.perm <- Y_data.perm$Outcom
          cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]
          set.seed(123)#set.seed(123)
          cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                                lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
          lamb <- cv_lasso$lambda.min
          opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                              family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
          bet.vec <- opt_model$beta[mediation_set]
          nde.perm[jk] <- opt_model$beta[(1+p)]

          print(jk)
        }
        ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

        nie = niea = niep = nie.p = niea.p = niep.p = ci.nie = ci.niea = ci.niep = NA}else{
          #######confidence interval ######
          nde.perm <- c()
          niea.perm <- matrix(NA,n.perm,length(mediation_set))
          niep.perm <- matrix(NA,n.perm,length(mediation_set))
          nie.perm <- matrix(NA,n.perm,length(mediation_set))
          Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)

          for(jk in 1:n.perm){
            set.seed(jk+2000)
            Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
            tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
            M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
            t.perm <- Y_data.perm$Exposure
            Y_mat.perm <- Y_data.perm$Outcom
            cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]

            set.seed(123)#set.seed(123)
            cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                                  lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
            lamb <- cv_lasso$lambda.min
            opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                                family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
            bet.vec <- opt_model$beta[mediation_set]
            nde.perm[jk] <- opt_model$beta[(1+p)]

            for(j in 1:length(mediation_set)){
              Y_data.dcrt <- as.data.frame(cbind(m.fi = as.vector(M_mat.perm[,mediation_set[j]]),Exposure=t.perm,Outcom=Y_mat.perm,Cov=cov.perm))

              zip.fit <-  try(zeroinfl(m.fi ~ Exposure+Cov | Exposure+Cov,offset=tij_mat.perm, data = Y_data.dcrt,dist = c("poisson"),link = c("logit")), silent = TRUE)
              if(inherits(zip.fit,"try-error")){niea.perm[jk,j]=niep.perm[jk,j]=nie.perm[jk,j]=NA}
              else{
                beta.val <- bet.vec[j]
                alpha.val<- summary(zip.fit)$coefficients$count[1:(2+ncol(cov)),1]
                gamma.val <- summary(zip.fit)$coefficients$zero[1:(2+ncol(cov)),1]


                niea.perm[jk,j] <- mean(beta.val*(1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(gamma.val[-2]))))*
                                          (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(alpha.val)+tij_mat.perm)-
                                             exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(alpha.val[-2])+tij_mat.perm)))
                niep.perm[jk,j] <- mean(beta.val*(exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(alpha.val)+tij_mat.perm))*
                                          ((1/(1+exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(gamma.val))))-
                                             (1/(1+exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(gamma.val[-2]))))))

                nie.perm[jk,j] <-niea.perm[jk,j]+niep.perm[jk,j]

              }

            }

            print(jk)
          }


          if(length(mediation_set)==1){
            ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.niea <- quantile(as.vector(niea.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.niep <- quantile(as.vector(niep.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.nie <- quantile(as.vector(nie.perm), probs = c(0.025, 0.975),na.rm=TRUE)
          }else{
            ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

            ci.niea <- apply(niea.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
            ci.niep <- apply(niep.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
            ci.nie <- apply(nie.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
          }

        }
      return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie, NIEA = niea, NIEP = niep,
                  NDE.pval = nde.p, NIE.pval = nie.p, NIEA.pval = niea.p, NIEP.pval = niep.p,
                  NDE.CI = ci.nde, NIE.CI = ci.nie, NIEA.CI = ci.niea, NIEP.CI = ci.niep,
                  AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
      ))
    }else{
      return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie, NIEA = niea, NIEP = niep,
                  NDE.pval = nde.p, NIE.pval = nie.p, NIEA.pval = niea.p, NIEP.pval = niep.p,
                  AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
      ))
    }
  }

  if(modelm=="NB"){
    for(j in 1:length(selection_set)){
      library(MASS)
      Y_data <- data.frame(Exposure=Exposure,Mediator=as.vector(M_mat[,selection_set[j]]),cov=cov)

        nb.fit <-  glm.nb(Mediator ~ Exposure+cov, data = Y_data)
        residualsm[[j]] <- nb.fit$residuals
        aic[j] <-  AIC(nb.fit); bic[j] <-  BIC(nb.fit)
        final.p[j] <- summary(nb.fit)$coefficients[2,4]
        alpha_mat[j,]<- summary(nb.fit)$coefficients[1:(ncol(cov)+2),1]
        alpha.p[j] <- summary(nb.fit)$coefficients[2,4]

        beta.val <- beta_fit[selection_set[j]];
        alpha.val <- alpha_mat[j,];

        nie[j] <- mean(beta.val*
                        (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val)+tij_mat)-
                           exp(as.matrix(data.frame(rep(1,nrow(M_mat)),cov))%*%matrix(alpha.val[-2])+tij_mat)))

}
    bind.p <- rbind(p.adjust(index.p,method),p.adjust(final.p,method))
    joint.p <- apply(bind.p, 2, max)
    index.mi <- which(joint.p <= FDR)
    mediation_set <- selection_set[index.mi]

    nie <- nie[index.mi]
    niea <- niea[index.mi]
    niep <- niep[index.mi]

    nie.p <- joint.p[index.mi]

    bind.p <- rbind(p.adjust(index.p,method),p.adjust(alpha.p,method))
    joint.p <- apply(bind.p, 2, max)
    niea.p <- joint.p[index.mi]

    bind.p <- rbind(p.adjust(index.p,method),p.adjust(gamma.p,method))
    joint.p <- apply(bind.p, 2, max)
    niep.p <- joint.p[index.mi]
    if(CI==TRUE){
      if(length(mediation_set)==0){
        nde.perm <- c()
        Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)
        for(jk in 1:n.perm){
          set.seed(jk+2000)
          Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
          tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
          M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
          t.perm <- Y_data.perm$Exposure
          Y_mat.perm <- Y_data.perm$Outcom
          cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]
          set.seed(123)#set.seed(123)
          cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                                lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
          lamb <- cv_lasso$lambda.min
          opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                              family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
          bet.vec <- opt_model$beta[mediation_set]
          nde.perm[jk] <- opt_model$beta[(1+p)]

          print(jk)
        }
        ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

        nie = niea = niep = nie.p = niea.p = niep.p = ci.nie = ci.niea = ci.niep = NA}else{
          #######confidence interval ######
          nde.perm <- c()
          niea.perm <- matrix(NA,n.perm,length(mediation_set))
          niep.perm <- matrix(NA,n.perm,length(mediation_set))
          nie.perm <- matrix(NA,n.perm,length(mediation_set))
          Y_data <- data.frame(M_mat,Cov=cov,Exposure=Exposure,Outcom=Y)

          for(jk in 1:n.perm){
            set.seed(jk+2000)
            Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]
            tij_mat.perm <-log(rowSums(Y_data.perm[,1:ncol(M_mat)]))
            M_mat.perm <- Y_data.perm[,1:ncol(M_mat)]
            t.perm <- Y_data.perm$Exposure
            Y_mat.perm <- Y_data.perm$Outcom
            cov.perm <- Y_data.perm[,(ncol(M_mat)+1):(ncol(M_mat)+ncol(cov))]

            set.seed(123)#set.seed(123)
            cv_lasso <- cv.glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, family = modely,
                                  lambda = NULL, dfmax = as.integer(ncol(M_mat.perm) / 2))
            lamb <- cv_lasso$lambda.min
            opt_model <- glmnet(as.matrix(data.frame(M_mat.perm,t.perm,cov.perm)), Y_mat.perm, alpha = 1, lambda = lamb,
                                family = modely, dfmax = as.integer(ncol(M_mat.perm) / 2))
            bet.vec <- opt_model$beta[mediation_set]
            nde.perm[jk] <- opt_model$beta[(1+p)]

            for(j in 1:length(mediation_set)){
              Y_data.dcrt <- as.data.frame(cbind(m.fi = as.vector(M_mat.perm[,mediation_set[j]]),Exposure=t.perm,Outcom=Y_mat.perm,Cov=cov.perm))

              nb.fit <-  try(glm.nb(m.fi ~ Exposure+Cov, data = Y_data.dcrt), silent = TRUE)
              if(inherits(nb.fit,"try-error")){niea.perm[jk,j]=niep.perm[jk,j]=nie.perm[jk,j]=NA}
              else{
                beta.val <- bet.vec[j]
                alpha.val<- summary(nb.fit)$coefficients$count[1:(2+ncol(cov)),1]


                nie.perm[jk,j] <- mean(beta.val*
                                          (exp(as.matrix(data.frame(1,rep(1,nrow(M_mat)),cov.perm))%*%matrix(alpha.val)+tij_mat.perm)-
                                             exp(as.matrix(data.frame(rep(1,nrow(M_mat.perm)),cov.perm))%*%matrix(alpha.val[-2])+tij_mat.perm)))



              }

            }

            print(jk)
          }


          if(length(mediation_set)==1){
            ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.niea <- quantile(as.vector(niea.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.niep <- quantile(as.vector(niep.perm), probs = c(0.025, 0.975),na.rm=TRUE)
            ci.nie <- quantile(as.vector(nie.perm), probs = c(0.025, 0.975),na.rm=TRUE)
          }else{
            ci.nde <- quantile(as.vector(nde.perm), probs = c(0.025, 0.975),na.rm=TRUE)

            ci.niea <- apply(niea.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
            ci.niep <- apply(niep.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
            ci.nie <- apply(nie.perm,2,function(x){
              quantile(x, probs = c(0.025, 0.975))
            })
          }

        }
      return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie,
                  NDE.pval = nde.p, NIE.pval = nie.p,
                  NDE.CI = ci.nde, NIE.CI = ci.nie,
                  AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
      ))
    }else{
      return(list(Mediators = mediation_set, NDE = beta_2, NIE = nie,
                  NDE.pval = nde.p, NIE.pval = nie.p,
                  AIC = aic, BIC = bic, residualsy=residualsy, residualsm=residualsm
      ))
    }
  }



  }
}



