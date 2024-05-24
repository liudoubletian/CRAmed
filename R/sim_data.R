#' Simulate the ZINB data
#'
#' @param otu_n the number od taxa.
#' @param num the number of samples.
#' @return Returns the data list, including the treatment, mediators and outcome.
#' @examples
#' otu_n <- 50;num <- 50
#' set.seed(123)
#' sim_zinb.mat <- sim_zinb(otu_n,num)
#' @export

sim_zinb <- function(otu_n,num){
  beta <- rep(0,otu_n+2)
  beta[1] <- 1
  beta[2] <- -2
  beta[3:12] <- rnorm(10,0.2,1)

  alphak <-rep(0,otu_n)
  gammak <- rep(0,otu_n)
  alpha.rand <- rnorm(25,-0.5,1)
  gamma.rand <- rnorm(25,-1,1)

  alphak[1:4] <- alpha.rand[1:4]
  alphak[11:15] <- alpha.rand[11:15]
  alphak[21:25] <- alpha.rand[21:25]

  gammak[1:3] <- gamma.rand[1:3]
  gammak[5] <- gamma.rand[5]
  gammak[11:20] <- gamma.rand[11:20]

  M_mat <- matrix(NA,num,otu_n)
  tij_mat <- matrix(NA,num,otu_n)

  t=matrix(rep(c(1,0),each=num/2),ncol=1)
  x=as.matrix(data.frame(itc=1,t))

  for(j in 1:otu_n){
  set.seed(1+j)
  ck <- runif(1,-2,2)
  logit.p  <- x %*% matrix(c(ck,gammak[j]))
  p  <- 1 / (1 + exp(-logit.p))
  ck <- -7
  tij_mat[,j] <- runif(num,7.1,10.5)
  lamda=exp(x%*%matrix(c(ck,alphak[j]))+tij_mat[,j])
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
Y_mat <- as.matrix(data.frame(x,M_mat))%*%matrix(beta)+rnorm(num)
return(list(M_mat=M_mat,Y_mat=Y_mat,t=t))
}
