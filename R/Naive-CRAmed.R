#' Naive.CRAmed: A conditional randomization test for normal microbiome data
#'
#' @param M_mat a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param Y a continuous outcome.
#' @param Exposure a binary exposure.
#' @param FDR the threshold value for p-values.
#' @param method the correction method for multiple testing.

#' @return Return the identified mediators and corresponding p-values of NIE.
#' @examples
#' otu_n <- 50;num <- 50
#' set.seed(12)
#' sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)
#' naive.cramed.res <- Naive.CRAmed(M_mat=log(sim_zinb.mat$M_mat+1), Y=sim_zinb.mat$Y, Exposure=sim_zinb.mat$trt)
#' @export Naive.CRAmed
#' @export Creat_condition_gaussian
#' @export dCRT
#' @importFrom glmnet  cv.glmnet
#' @import MASS
#' @import plyr
#' @import pscl

Creat_condition_gaussian <- function(X, indx, Sigma = NULL, lambda = NULL){
  library(glmnet)
  if (is.null(Sigma)){
    p <- length(X[1,])
    library(glmnet)
    cv_lasso_x <- cv.glmnet(X[,-indx], X[,indx], alpha = 1,
                            intercept = T, dfmax = as.integer(p / 2))
    print(3)
    lambda_x <- cv_lasso_x$lambda.min
    opt_model_x <- glmnet(X[,-indx], X[,indx], alpha = 1,
                          lambda = lambda_x, intercept = T, dfmax = as.integer(p / 2))
    beta_x <- opt_model_x$beta
    X_bar <- predict(opt_model_x, X[,-indx])
    x_res <- X[,indx] - X_bar
    sigma2_x <- mean(x_res^2)
  }else{
    beta_x <- solve(Sigma[-indx, -indx], Sigma[-indx, indx])
    X_bar <- X[ ,-indx] %*% beta_x
    sigma2_x <- Sigma[indx, indx] - Sigma[indx, -indx] %*% beta_x
    ##print(beta_x)
  }
  return(list(mean_x = X_bar, sigma2_x = sigma2_x, gamma = beta_x))
}


dCRT <-function(X, Y, Sigma_X = NULL, FDR = 0.05, candidate_set ,
                mean_X_mat = NULL, Gen_X = NULL, CDR = 'Consist w model',
                model = 'Gaussian_lasso', d.interaction = F, RF.num.trees = c(100, 30),
                lambda.seq = NULL, k = NULL, M = 50000,
                central = F, MC_free = F,  X_quantile_mat = NULL,
                X_cond_var_mat = NULL, return.pvalues = T,
                eps_dI = 1e-8) {


  library(glmnet)
  library("randomForest")
  if (!model %in% c('Gaussian_lasso', 'Binomial_lasso', 'RF')){
    print('Error: please correctly specify the model for Y.')
    return()
  }

  if (model == 'Gaussian_lasso'){
    model <- 'gaussian'
  }
  if (model == 'Binomial_lasso'){
    model <- 'binomial'
  }

  if (is.null(Gen_X) & ! is.null(mean_X_mat)){
    print('Error: please specify distribution of X or centralize the gaussian X.')
    return()
  }

  if (! is.null(Gen_X) & is.null(mean_X_mat)){
    print('Error: please specify the conditional mean of X')
    return()
  }
  if (MC_free == T){
    type = 'cov'
  }else{
    type = 'beta'
  }

  p <- length(X[1,])
  n <- length(Y)
  if (is.null(k)){
    if (d.interaction == T){
      k = min(as.integer(2 * log(p)), p - 1)
    }else{
      k = 0
    }
  }

  if (central == T){
    mean_X_mat <- mean_X_mat - colMeans(X)
    X <- X - colMeans(X)
  }

  ############### CDR ###############

  if (CDR == 'Consist w model'){
    if (model == 'gaussian' | model == 'binomial'){
      cv_lasso <- cv.glmnet(X, Y, alpha = 1, family = model,
                            lambda = lambda.seq, dfmax = as.integer(p / 2))
      lamb <- cv_lasso$lambda.min
      opt_model <- glmnet(X, Y, alpha = 1, lambda = lamb,
                          family = model, dfmax = as.integer(p / 2))
      beta_fit <- as.numeric(opt_model$beta)
      selection_set <- which(beta_fit != 0)
    }
  }else{
    if (CDR == 'No'){
      selection_set <- 1:p
    }else{
      selection_set <- CDR(Y, X)
    }
  }

  selection_set <- intersect(candidate_set, selection_set)
  print('############# Candidate set ###############')
  print(selection_set)


  ############## dCRT ##############

  sigma2_lst <- vector('list', length(selection_set))
  X_bar_lst <- vector('list', length(selection_set))
  offsets_lst <- vector('list', length(selection_set))
  imp_obe_lst <- vector('list', length(selection_set))
  X_use_lst <- vector('list', length(selection_set))
  pvl_lst <- rep(1, p)
  est_lst <- c()
  print(1)
  for (j in 1:length(selection_set)){

    indx <- selection_set[j]

    ############### Distill X ##################
    if (is.null(Gen_X)){
      library(glmnet)
      Cond_X <- Creat_condition_gaussian(X, indx, Sigma = Sigma_X)
      print(2)
      mean_X <- Cond_X$mean_x
      sigma2_x <- Cond_X$sigma2_x
      sigma2_lst[[j]] <- sigma2_x
      X_bar_lst[[j]] <- mean_X
    }else{
      mean_X <- mean_X_mat[,indx]
      X_bar_lst[[j]] <- mean_X
    }

    X_res_ob <- X[,indx] - mean_X

    ############### Distill Y and Importance ##################

    if (model == 'RF'){
      n_tree <- RF.num.trees[1]
      n_tree_imp <- RF.num.trees[2]

      if (k >= 1){
        rf.fit <- randomForest(x = X[,-indx], y = Y, ntree = n_tree)
        offsets <- rf.fit$predicted
        offsets_lst[[j]] <- offsets

        imp_fit <- as.vector(rf.fit$importance)
        imp_sort <- sort(imp_fit, decreasing = T, index.return = T)
        index_use <- imp_sort$ix[1:k]

        X_use <- X[,-indx][,index_use]
        X_use_lst[[j]] <- X_use

        # Observed importance

        rf.indx.fit <- randomForest(x = as.matrix(cbind(X[,indx] - mean_X, offsets, X_use)),
                                    Y, ntree = n_tree_imp)
        imp_obe <- rf.indx.fit$importance[1]
        imp_obe_lst[[j]] <- imp_obe

      }else{
        rf.fit <- randomForest(x = X[,-indx] - c(mean_X), y = Y, ntree = n_tree)
        eps_res <- Y - rf.fit$predicted

        # Observed importance
        if (MC_free == F | is.null(Gen_X)){
          imp_obe <- abs(mean(X_res_ob * eps_res))
        }
      }

    }

    if (model == 'binomial' | model == 'gaussian'){
      if (k >= 1){
        if (model == 'gaussian'){
          family = gaussian()
        }
        if (model == 'binomial'){
          family = binomial()
        }
        cv_lasso <- cv.glmnet(X[,-indx], Y, alpha = 1, family = model,
                              lambda = lambda.seq, dfmax = as.integer(p / 2))
        lamb <- cv_lasso$lambda.min
        opt_model <- glmnet(X[,-indx], Y, alpha = 1, lambda = lamb,
                            family = model, dfmax = as.integer(p / 2))

        beta_leave <- opt_model$beta
        beta_sort <- sort(abs(as.vector(beta_leave)), decreasing = T, index.return = T)
        index_use <- beta_sort$ix[1:k]
        X_use <- X[,-indx][,index_use]
        X_use_lst[[j]] <- X_use

        offsets <- as.vector(opt_model$a0 + X[,-indx] %*% opt_model$beta)

        if (model == 'binomial'){
          eps_res <- Y - 1 / (1 + exp(- offsets))
        }
        if (model == 'gaussian'){
          eps_res <- Y - offsets
        }

        offsets_lst[[j]] <- offsets

        print('############# distilled ###############')
        print(indx)

        if (type == 'cov'){

          if (!is.null(Gen_X)){
            X_quant <- X_quantile_mat[,indx]
            X_res_ob <- qnorm(X_quant, mean = rep(0, n), sd = 1)
            sigma2_x <- 1
          }

          weight_inter <- 1 / sqrt(k)
          W_inter <- eps_res
          for (l in 1:length(X_use[1,])){
            W_inter <- cbind(W_inter, X_use[,l] * eps_res)
          }
          XTX <- t(cbind(1, X_use)) %*% cbind(1, X_use)
          XTX_inv <- solve(XTX)

          Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% X_res_ob
          WTW <- t(W_inter) %*% W_inter
          svd_WTW <- svd(WTW)
          root_WTW <- svd_WTW$u %*% diag(sqrt(svd_WTW$d)) %*% t(svd_WTW$v)
          lambda_W <- svd(root_WTW %*% XTX_inv %*%
                            diag(c(1, rep(weight_inter^2, k))) %*% XTX_inv %*% root_WTW)$d
          imp_obe_lst[[j]] <- sum(Z_dI^2) / c(sigma2_x)
        }

        if (type == 'beta'){
          weight_inter <- 1 / sqrt(k)
          W_inter <- cbind(1, X_use) * as.vector(eps_res)
          W_inter_X <- cbind(1, X_use) * as.vector(X_res_ob)
          XTX <- t(W_inter_X) %*% W_inter_X
          XTX_inv <- solve(XTX)
          Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% X_res_ob
          imp_obe_lst[[j]] <- sum(Z_dI^2)

        }

      }else{
        cv_lasso_null <- cv.glmnet(X[,-indx], Y, alpha = 1, family = model,
                                   lambda = lambda.seq, dfmax = as.integer(p / 2))
        lamb_null <- cv_lasso_null$lambda.min
        model_res_null <- glmnet(X[,-indx], Y, alpha = 1, lambda = lamb_null,
                                 family = model, dfmax = as.integer(p / 2))
        if (model == 'binomial'){
          eps_res <- Y - 1 / (1 + exp(- predict(model_res_null, X[,-indx])))
        }
        if (model == 'gaussian'){
          eps_res <- Y - predict(model_res_null, X[,-indx])
        }
        if (MC_free == F | is.null(Gen_X)){
          if (type == 'cov'){
            imp_obe <- abs(mean(X_res_ob * eps_res))
          }
          if (type == 'beta'){
            imp_obe <- abs(mean(X_res_ob * eps_res)) / mean(X_res_ob^2)
          }
        }
        print('############# distilled ###############')
        print(indx)
      }
    }


    ############### Estimate p-values ##################

    if (k == 0){
      if (is.null(Gen_X)){

        if (type == 'cov'){
          emp_var <- mean(eps_res^2) * sigma2_x
          pvl <- 2 * pnorm(- sqrt(n) * abs(imp_obe) / sqrt(emp_var))
          pvl_lst[indx] <- pvl
        }
        if (type == 'beta'){
          X_res_sample <- rnorm(n * M, 0, sqrt(sigma2_x))
          X_res_sample <- matrix(X_res_sample, n, M)
          var_lst_sample <- unlist(lapply(c(1:M), function(j){
            mean((X_res_sample[,j])^2) }))
          t_lst <- abs(t(X_res_sample) %*% eps_res / n) / var_lst_sample
          pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
        }

      }
      else{
        if (MC_free == T){
          if (is.null(X_quantile_mat) | is.null(X_cond_var_mat)){
            print('Please specify conditional quantile of X when using MCf-d0CRT')
            return()
          }
          X_quant <- X_quantile_mat[,indx]
          X_cond_var <- X_cond_var_mat[,indx]
          X_norm <- qnorm(X_quant, mean = rep(0, n), sd = sqrt(X_cond_var))
          imp_obe <- abs(mean(X_norm * eps_res))
          emp_var <- mean(X_cond_var * eps_res^2)
          pvl <- 2 * pnorm(- sqrt(n) * abs(imp_obe) / sqrt(emp_var))
          pvl_lst[indx] <- pvl

        }else{
          X_resample <- Gen_X(X, indx, num = M)
          X_res_resample <- X_resample - mean_X

          if (type == 'cov'){
            t_lst <- t(X_res_resample) %*% eps_res / n
            t_lst <- c(imp_obe, abs(t_lst))
          }
          if (type == 'beta'){
            t_lst <- unlist(lapply(c(1:M), function(j){
              mean(X_res_resample[,j] * eps_res) / mean((X_res_resample[,j])^2)
            }))
            t_lst <- c(imp_obe, abs(t_lst))
          }
          pvl_lst[indx] <- mean(ifelse(t_lst >= imp_obe, 1, 0))
        }
      }

      print('############# d0CRT pvalue ###############')
      print(paste(indx, pvl_lst[indx], sep = ': '))

    }else{

      if (model == 'RF'){
        Z_vec <- c(1)
        for (m in 1:M) {
          if (is.null(Gen_X)){
            delta_gen <- rnorm(n, 0, sqrt(sigma2_x))
            X_sample <- mean_X + delta_gen
          }else{
            X_sample <- Gen_X(X, indx, num = 1)
          }
          rf.sam.fit <- randomForest(x = as.matrix(cbind(X_sample - mean_X, offsets, X_use)),
                                     Y, ntree = n_tree_imp)
          imp_resample <- rf.sam.fit$importance[1]
          Z <- ifelse(imp_resample >= imp_obe, 1, 0)
          Z_vec <- c(Z_vec, Z)
        }
        pvl_lst[indx] <- mean(Z_vec)
      }else{
        if (type == 'cov'){
          pvl <- imhof(imp_obe_lst[[j]], lambda_W, epsabs = eps_dI)
          pvl_lst[indx] <- abs(pvl$Qq)
        }

        if (type == 'beta'){

          if (is.null(Gen_X)){
            X_res_sample <- rnorm(n * M, 0, sqrt(sigma2_x))
            X_res_sample <- matrix(X_res_sample, n, M)
          }else{
            X_resample <- Gen_X(X, indx, num = M)
            X_res_sample <- X_resample - mean_X
          }

          W_inter <- cbind(1, X_use) * as.vector(eps_res)
          weight_inter <- 1 / sqrt(k)

          t_lst <- unlist(lapply(c(1:M), function(j){
            W_inter_X <- cbind(1, X_use) * as.vector(X_res_sample[,j])
            XTX <- t(W_inter_X) %*% W_inter_X
            XTX_inv <- solve(XTX)
            Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% X_res_sample[,j]
            sum(Z_dI^2)
          }))

          pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe_lst[[j]], 1, 0)))

        }
      }

      print('############# dkCRT pvalue ###############')
      print(paste(indx, pvl_lst[indx], sep = ': '))

    }
  }

  CRT_BHq_lst <- pvl_lst
  selection_set_CRT <- which(CRT_BHq_lst <= FDR)
  if(model == 'RF'){
    if (return.pvalues == T){
      return(list(select_set = selection_set_CRT, p_values = pvl_lst[selection_set_CRT],p.total=pvl_lst))
    }else{
      return(list(select_set = selection_set_CRT))
    }
  }
  else{
    if (return.pvalues == T){
      return(list(select_set = selection_set_CRT, p_values = pvl_lst[selection_set_CRT],est = beta_fit[selection_set_CRT],est.total = beta_fit,lasso.sel = selection_set,p.total=pvl_lst))
    }else{
      return(list(select_set = selection_set_CRT))
    }
  }
}
Naive.CRAmed <- function(M_mat, Y, Exposure, FDR = 0.05, method="fdr") {



  dcrt.res <- dCRT(X=cbind(log(as.matrix(M_mat)+1),Exposure), Y=Y, candidate_set = 1:ncol(M_mat))


  ID_test <- dcrt.res$lasso.sel
  index.p <- dcrt.res$p.total

  mat.alpha <-matrix(NA,length(ID_test),2)

  final.p <- c()
  for(j in 1:length(ID_test)){
    Y_data <- as.data.frame(cbind(Exposure,log(M_mat[,ID_test[j]]+1),Y))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    final.p[j] <-  summary(lm(Mediator ~ Treat, data = Y_data))$coefficients[2,4]
  }

  bind.p <- rbind(p.adjust(index.p,method),p.adjust(final.p,method))
  joint.p <- apply(bind.p, 2, max)

  mediation_set <- ID_test[which(joint.p<FDR)]
  nie.pval <- joint.p[which(joint.p<FDR)]
  return(list(Mediators = mediation_set, NIE.pval = nie.pval ))
  }

