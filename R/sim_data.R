#' Simulate the ZINB data
#'
#' @param otu_n the number od taxa.
#' @param num the number of samples.
#' @param alpha the corresponding regression coefficients for alpha.
#' @param beta the corresponding regression coefficients for beta.
#' @param gamma the corresponding regression coefficients for gamma.
#' @param conf.num the number of confounders.

#' @return Returns the data list, including the treatment, mediators and outcome.
#' @examples
#' otu_n <- 100;num <- 50
#' set.seed(123)
#' sim_zinb.mat <- sim_zinb(otu_n, num, alpha=-2, beta=2, gamma=-2)
#' @import aster
#' @export

sim_zinb <- function(otu_n,num,alpha,beta,gamma,conf.num=NULL){
  beta.vec <- rep(0,otu_n+2)
  beta.vec[1] <- 1
  beta.vec[2] <- -2
  beta.vec[3:12] <- rnorm(10,beta,1)

  alpha.vec <-rep(0,otu_n)
  gamma.vec <- rep(0,otu_n)
  alpha.rand <- rnorm(25,alpha,1)
  gamma.rand <- rnorm(25,gamma,1)

  alpha.vec[1:4] <- alpha.rand[1:4]
  alpha.vec[11:15] <- alpha.rand[11:15]
  alpha.vec[21:25] <- alpha.rand[21:25]

  gamma.vec[1:3] <- gamma.rand[1:3]
  gamma.vec[5] <- gamma.rand[5]
  gamma.vec[11:20] <- gamma.rand[11:20]

  M_mat <- matrix(NA,num,otu_n)
  tij_mat <- matrix(NA,num,otu_n)

  trt=matrix(rep(c(1,0),each=num/2),ncol=1)
  x=as.matrix(data.frame(itc=1,trt))

  if(is.null(conf.num)){
  for(j in 1:otu_n){
  set.seed(1+j)
  intercept <- runif(1,-2,2)
  logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))
  p  <- 1 / (1 + exp(-logit.p))
  intercept <- -7
  tij_mat[,j] <- runif(num,7.1,10.5)
  lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j])
  M1 <- rbinom(num, 1, p)
  for (i in 1:num) {
    if (M1[i] == 1) {
      M_mat[i,j] = 0
    }
    if (M1[i] == 0) {
      phi <- runif(1,0.1,10)
      M_mat[i,j] = rnbinom(1, mu=lamda[i],size=phi)
    }
  }
}

set.seed(1)
Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+rnorm(num)
  }else{
    X_con <- matrix(rnorm(conf.num*num),ncol=conf.num)
    alpha.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    gamma.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    for(j in 1:otu_n){
      set.seed(1+j)
      intercept <-runif(1,-2,2)
      logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))+X_con%*%gamma.con[,j]
      p  <- 1 / (1 + exp(-logit.p))
      intercept <- -7
      tij_mat[,j] <- runif(num,7.1,10.5)
      lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j]+X_con%*%alpha.con[,j])
      M1 <- rbinom(num, 1, p)
      for (i in 1:num) {
        if (M1[i] == 1) {
          M_mat[i,j] = 0
        }
        if (M1[i] == 0) {
          phi <- runif(1,0.1,10)
          M_mat[i,j] = rnbinom(1, mu=lamda[i],size=phi)
        }
      }
    }
    set.seed(1)
    beta.con <- matrix(rnorm(conf.num))
    Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+X_con%*%beta.con+rnorm(num)
  }
return(list(M_mat=M_mat,Y=Y,trt=trt))
}
sim_hnb <- function(otu_n,num,alpha,beta,gamma,conf.num=NULL){
  library(aster)
  beta.vec <- rep(0,otu_n+2)
  beta.vec[1] <- 1
  beta.vec[2] <- -2
  beta.vec[3:12] <- rnorm(10,beta,1)

  alpha.vec <-rep(0,otu_n)
  gamma.vec <- rep(0,otu_n)
  alpha.rand <- rnorm(25,alpha,1)
  gamma.rand <- rnorm(25,gamma,1)

  alpha.vec[1:4] <- alpha.rand[1:4]
  alpha.vec[11:15] <- alpha.rand[11:15]
  alpha.vec[21:25] <- alpha.rand[21:25]

  gamma.vec[1:3] <- gamma.rand[1:3]
  gamma.vec[5] <- gamma.rand[5]
  gamma.vec[11:20] <- gamma.rand[11:20]

  M_mat <- matrix(NA,num,otu_n)
  tij_mat <- matrix(NA,num,otu_n)

  trt=matrix(rep(c(1,0),each=num/2),ncol=1)
  x=as.matrix(data.frame(itc=1,trt))

  if(is.null(conf.num)){
    for(j in 1:otu_n){
      set.seed(1+j)
      intercept <- runif(1,-2,2)
      logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))
      p  <- 1 / (1 + exp(-logit.p))
      intercept <- -7
      tij_mat[,j] <- runif(num,7.1,10.5)
      lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j])
      M1 <- rbinom(num, 1, p)
      for (i in 1:num) {
        if (M1[i] == 1) {
          M_mat[i,j] = 0
        }
        if (M1[i] == 0) {
          phi <- runif(1,0.1,10)
          M_mat[i,j] =  rktnb(n=1, size=phi,k=0 , mu=lamda[i], xpred=1)
        }
      }
    }

    set.seed(1)
    Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+rnorm(num)
  }else{
    X_con <- matrix(rnorm(conf.num*num),ncol=conf.num)
    alpha.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    gamma.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    for(j in 1:otu_n){
      set.seed(1+j)
      intercept <-runif(1,-2,2)
      logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))+X_con%*%gamma.con[,j]
      p  <- 1 / (1 + exp(-logit.p))
      intercept <- -7
      tij_mat[,j] <- runif(num,7.1,10.5)
      lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j]+X_con%*%alpha.con[,j])
      M1 <- rbinom(num, 1, p)
      for (i in 1:num) {
        if (M1[i] == 1) {
          M_mat[i,j] = 0
        }
        if (M1[i] == 0) {
          phi <- runif(1,0.1,10)
          M_mat[i,j] = rnbinom(1, mu=lamda[i],size=phi)
        }
      }
    }
    set.seed(1)
    beta.con <- matrix(rnorm(conf.num))
    Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+X_con%*%beta.con+rnorm(num)
  }
  return(list(M_mat=M_mat,Y=Y,trt=trt))
}


sim_hp <- function(otu_n,num,alpha,beta,gamma,conf.num=NULL){
  library(aster)
  beta.vec <- rep(0,otu_n+2)
  beta.vec[1] <- 1
  beta.vec[2] <- -2
  beta.vec[3:12] <- rnorm(10,beta,1)

  alpha.vec <-rep(0,otu_n)
  gamma.vec <- rep(0,otu_n)
  alpha.rand <- rnorm(25,alpha,1)
  gamma.rand <- rnorm(25,gamma,1)

  alpha.vec[1:4] <- alpha.rand[1:4]
  alpha.vec[11:15] <- alpha.rand[11:15]
  alpha.vec[21:25] <- alpha.rand[21:25]

  gamma.vec[1:3] <- gamma.rand[1:3]
  gamma.vec[5] <- gamma.rand[5]
  gamma.vec[11:20] <- gamma.rand[11:20]

  M_mat <- matrix(NA,num,otu_n)
  tij_mat <- matrix(NA,num,otu_n)

  trt=matrix(rep(c(1,0),each=num/2),ncol=1)
  x=as.matrix(data.frame(itc=1,trt))

  if(is.null(conf.num)){
    for(j in 1:otu_n){
      set.seed(1+j)
      intercept <- runif(1,-2,2)
      logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))
      p  <- 1 / (1 + exp(-logit.p))
      intercept <- -7
      tij_mat[,j] <- runif(num,7.1,10.5)
      lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j])
      M1 <- rbinom(num, 1, p)
      for (i in 1:num) {
        if (M1[i] == 1) {
          M_mat[i,j] = 0
        }
        if (M1[i] == 0) {
          phi <- runif(1,0.1,10)
          M_mat[i,j] = rnbinom(1, mu=lamda[i],size=phi)
        }
      }
    }

    set.seed(1)
    Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+rnorm(num)
  }else{
    X_con <- matrix(rnorm(conf.num*num),ncol=conf.num)
    alpha.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    gamma.con <- matrix(rnorm(conf.num*otu_n),nrow=conf.num)
    for(j in 1:otu_n){
      set.seed(1+j)
      intercept <-runif(1,-2,2)
      logit.p  <- x %*% matrix(c(intercept,gamma.vec[j]))+X_con%*%gamma.con[,j]
      p  <- 1 / (1 + exp(-logit.p))
      intercept <- -7
      tij_mat[,j] <- runif(num,7.1,10.5)
      lamda=exp(x%*%matrix(c(intercept,alpha.vec[j]))+tij_mat[,j]+X_con%*%alpha.con[,j])
      M1 <- rbinom(num, 1, p)
      for (i in 1:num) {
        if (M1[i] == 1) {
          M_mat[i,j] = 0
        }
        if (M1[i] == 0) {
          M_mat[i,j] = rktp(n=1, k=0 , mu=lamda[i], xpred=1)
        }
      }
    }
    set.seed(1)
    beta.con <- matrix(rnorm(conf.num))
    Y <- as.matrix(data.frame(x,M_mat))%*%matrix(beta.vec)+X_con%*%beta.con+rnorm(num)
  }
  return(list(M_mat=M_mat,Y=Y,trt=trt))
}
