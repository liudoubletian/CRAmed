####load packages
library(aster)
library(HDMT)
library(mediation)
library(pscl);packageVersion("pscl")# ‘1.5.5’
library(ncvreg)
library(doParallel)
library(foreach)
library(iterators)
library(parallel)
library(locfdr)
library(MAZE)
library(plyr)
library(MGLM)
library(fdrtool)
library(glmnet)
library(mvtnorm)
library(hdi)
library(MultiMed)
library(HIMA)
library(permute)
library(vegan)
######
dif <- 1:5
cramed.rec=cramed.pre=c()
ikt.rec=ikt.pre=hima.rec=hima.pre=c()
hima2.rec=hima2.pre=hdma.rec=hdma.pre=c()
multm.rec=multm.pre=ldm.rec=ldm.pre=c()
for(cc in 1:100){
  #######generate simulated data
  generate_zinb <- sim_zinb(otu_n=1000, num=100, alpha=-5, beta=5, gamma=-5, conf.num=3)
  M_mat <- generate_zinb$M_mat
  Y <-generate_zinb$Y
  trt <-generate_zinb$trt
  #######CRAmed
  results.cra <- try(CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05, n.times=100,
                            method="fdr",CI=FALSE),silent = TRUE)

  dec.cra <-  results.cra$Mediators
  tp <- length(dec.cra[dec.cra%in%dif])
  fp <- length(dec.cra)-tp
  fn <- length(dif)-tp
  cramed.rec[cc] <- tp/length(dif)
  if(length(dec.cra)==0){
    cramed.pre[cc] =1
  }
  else{
    cramed.pre[cc] <- tp/(length(dec.cra))
  }
  #####IKT
  med.mat = matrix(NA,otu_n,2)
  for(j in 1:otu_n){
    Y_data <- as.data.frame(cbind(trt,log(M_mat[,j]+1),Y))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    b <-  lm(Mediator ~ Treat, data = Y_data)
    c <- lm(Outcom ~ Treat + Mediator, data=Y_data)
    ikt.res <- try(mediate(b, c, sims=100, treat="Treat", mediator="Mediator"), silent = TRUE)
    if(inherits(contcont, "try-error")) {med.mat[j,]=NA}
    else{
      ikt.res <- summary(ikt.res)
      med.mat[j,] <- c(ikt.res$d0.p,ikt.res$z0.p)
    }
    print(j)
  }

  ikt.cla <- which(p.adjust(med.mat[,1],"fdr")<0.05)
  tp <- length(ikt.cla[ikt.cla%in%dif])
  fp <- length(ikt.cla)-tp
  fn <- length(dif)-tp
  ikt.rec[cc] <- tp/length(dif)
  if(length(ikt.cla)==0 ){
    ikt.pre[cc] =1
  }
  else{
    ikt.pre[cc] <- tp/(length(ikt.cla))
  }

  #####HIMA
  M_mat <- as.data.frame(M_mat)
  colnames(M_mat) = seq(1:otu_n)
  hima.fit <- try(hima(X =trt, Y = Y, M =log(M_mat+1), COV.XM =NULL,scale = FALSE,M.family = c( "gaussian"),Y.family = c("gaussian"),
                       verbose = FALSE),silent=TRUE)
  if(inherits(hima.fit,"try-error")){hima.rec[cc]=hima.pre[cc]=NA}else{
    hima.sel <- hima.fit[which(hima.fit$BH.FDR<0.05),]
    hima.cla <- unlist(lapply(strsplit(rownames(hima.sel),'`'),function(x){x[2]}))

    tp <- length(hima.cla[hima.cla%in%dif])
    fp <- length(hima.cla)-tp
    fn <- length(dif)-tp
    hima.rec[cc] <- tp/length(dif)
    if(length(hima.cla)==0 ){
      hima.pre[cc] =1
    }
    else{
      hima.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HIMA2

  source("./HIMA2.R")
  hima2.fit <- try(HIMA2(X=as.matrix(trt),Y=Y,M=log(M_mat+1),Z=NULL,cut=0.05),silent = TRUE)
  if(inherits(hima2.fit,"try-error")){hima2.rec[cc]=hima2.pre[cc]=NA}
  else{
    hima2.cla <- unlist(lapply(strsplit(hima2.fit$M,'M'),function(x){x[2]}))
    tp <- length(hima2.cla[hima2.cla%in%dif])
    fp <- length(hima2.cla)-tp
    fn <- length(dif)-tp
    hima2.rec[cc] <- tp/length(dif)
    if(length(hima2.cla)==0 ){
      hima2.pre[cc] =1
    }
    else{
      hima2.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HDMA
  source("./utils.R")
  source("./hdma.R")

  hdma.res <-  try(hdma(trt, Y=Y, M=log(M_mat+1), COV.XM = NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN =NULL,
                        parallel = FALSE, ncore = 1, verbose = FALSE),silent=TRUE)
  if(inherits(hdma.res,"try-error")){
    hdma.rec[cc]=hdma.pre[cc]=NA
  }
  else{
    hdma.cla <- rownames(hdma.res)[which(p.adjust(hdma.res$P.value,"fdr")<0.05)]

    tp <- length(hdma.cla[hdma.cla%in%dif])
    fp <- length(hdma.cla)-tp
    fn <- length(dif)-tp
    hdma.rec[cc] <- tp/length(dif)
    if(length(hdma.cla)==0 ){
      hdma.pre[cc] =1
    }
    else{
      hdma.pre[cc] <- tp/(length(hdma.cla))
    }
  }

  #####Multimed
  set.seed(18)
  multm <- try(medTest(E=as.matrix(t),M=as.matrix(log(M_mat+1)),Y=Y_mat,Z=NULL,  w=1),silent = TRUE)
  if(inherits(multm,"try-error")){multm.rec[cc]=multm.pre[cc]=NA}else{
    multm.cla <- which(p.adjust(multm[,2],"fdr") <0.05)
    tp <- length(multm.cla[multm.cla%in%dif])
    fp <- length(multm.cla)-tp
    fn <- length(dif)-tp
    multm.rec[cc] <- tp/length(dif)
    if(length(multm.cla)==0 ){
      multm.pre[cc] =1
    }
    else{
      multm.pre[cc] <- tp/(length(multm.cla))
    }
  }







  #####LDM
  source("./LDM_fun.R")
  meta <-data.frame(trt=trt,Y=Y)
  ldm.res <- ldm(M_mat~ trt+Y,
                 data=meta, seed=67817,fdr.nominal=0.05,
                 test.mediation=TRUE)

  ldm.cla <- ldm.res$med.detected.otu.omni
  tp <- length(ldm.cla[ldm.cla%in%dif])
  fp <- length(ldm.cla)-tp
  fn <- length(dif)-tp
  ldm.rec[cc] <- tp/length(dif)
  if(length(ldm.cla)==0){
    ldm.pre[cc] =1
  }
  else{
    ldm.pre[cc] <- tp/(length(dec.cla))
  }





}

write.table(rbind(cramed.rec,cramed.pre,
                  ikt.rec,ikt.pre,hima.rec,hima.pre,
                  hima2.rec,hima2.pre,hdma.rec,hdma.pre,
                  multm.rec,multm.pre,
                  ldm.rec,ldm.pre), "./figure3_a.txt")


#######
dif <- 1:5
cramed.rec=cramed.pre=c()
ikt.rec=ikt.pre=hima.rec=hima.pre=c()
hima2.rec=hima2.pre=hdma.rec=hdma.pre=c()
multm.rec=multm.pre=ldm.rec=ldm.pre=c()
for(cc in 1:100){
  #######generate simulated data
  generate_hp <- sim_hp(otu_n=1000, num=200, alpha=-5, beta=5, gamma=-5, conf.num=3)
  M_mat <- generate_hp$M_mat
  Y <-generate_hp$Y
  trt <-generate_hp$trt
  #######CRAmed
  results.cra <- try(CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05, n.times=100,
                            method="fdr",CI=FALSE),silent = TRUE)

  dec.cra <-  results.cra$Mediators
  tp <- length(dec.cra[dec.cra%in%dif])
  fp <- length(dec.cra)-tp
  fn <- length(dif)-tp
  cramed.rec[cc] <- tp/length(dif)
  if(length(dec.cra)==0){
    cramed.pre[cc] =1
  }
  else{
    cramed.pre[cc] <- tp/(length(dec.cra))
  }
  #####IKT
  med.mat = matrix(NA,otu_n,2)
  for(j in 1:otu_n){
    Y_data <- as.data.frame(cbind(trt,log(M_mat[,j]+1),Y))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    b <-  lm(Mediator ~ Treat, data = Y_data)
    c <- lm(Outcom ~ Treat + Mediator, data=Y_data)
    ikt.res <- try(mediate(b, c, sims=100, treat="Treat", mediator="Mediator"), silent = TRUE)
    if(inherits(contcont, "try-error")) {med.mat[j,]=NA}
    else{
      ikt.res <- summary(ikt.res)
      med.mat[j,] <- c(ikt.res$d0.p,ikt.res$z0.p)
    }
    print(j)
  }

  ikt.cla <- which(p.adjust(med.mat[,1],"fdr")<0.05)
  tp <- length(ikt.cla[ikt.cla%in%dif])
  fp <- length(ikt.cla)-tp
  fn <- length(dif)-tp
  ikt.rec[cc] <- tp/length(dif)
  if(length(ikt.cla)==0 ){
    ikt.pre[cc] =1
  }
  else{
    ikt.pre[cc] <- tp/(length(ikt.cla))
  }

  #####HIMA
  M_mat <- as.data.frame(M_mat)
  colnames(M_mat) = seq(1:otu_n)
  hima.fit <- try(hima(X =trt, Y = Y, M =log(M_mat+1), COV.XM =NULL,scale = FALSE,M.family = c( "gaussian"),Y.family = c("gaussian"),
                       verbose = FALSE),silent=TRUE)
  if(inherits(hima.fit,"try-error")){hima.rec[cc]=hima.pre[cc]=NA}else{
    hima.sel <- hima.fit[which(hima.fit$BH.FDR<0.05),]
    hima.cla <- unlist(lapply(strsplit(rownames(hima.sel),'`'),function(x){x[2]}))

    tp <- length(hima.cla[hima.cla%in%dif])
    fp <- length(hima.cla)-tp
    fn <- length(dif)-tp
    hima.rec[cc] <- tp/length(dif)
    if(length(hima.cla)==0 ){
      hima.pre[cc] =1
    }
    else{
      hima.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HIMA2

  source("./HIMA2.R")
  hima2.fit <- try(HIMA2(X=as.matrix(trt),Y=Y,M=log(M_mat+1),Z=NULL,cut=0.05),silent = TRUE)
  if(inherits(hima2.fit,"try-error")){hima2.rec[cc]=hima2.pre[cc]=NA}
  else{
    hima2.cla <- unlist(lapply(strsplit(hima2.fit$M,'M'),function(x){x[2]}))
    tp <- length(hima2.cla[hima2.cla%in%dif])
    fp <- length(hima2.cla)-tp
    fn <- length(dif)-tp
    hima2.rec[cc] <- tp/length(dif)
    if(length(hima2.cla)==0 ){
      hima2.pre[cc] =1
    }
    else{
      hima2.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HDMA
  source("./utils.R")
  source("./hdma.R")

  hdma.res <-  try(hdma(trt, Y=Y, M=log(M_mat+1), COV.XM = NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN =NULL,
                        parallel = FALSE, ncore = 1, verbose = FALSE),silent=TRUE)
  if(inherits(hdma.res,"try-error")){
    hdma.rec[cc]=hdma.pre[cc]=NA
  }
  else{
    hdma.cla <- rownames(hdma.res)[which(p.adjust(hdma.res$P.value,"fdr")<0.05)]

    tp <- length(hdma.cla[hdma.cla%in%dif])
    fp <- length(hdma.cla)-tp
    fn <- length(dif)-tp
    hdma.rec[cc] <- tp/length(dif)
    if(length(hdma.cla)==0 ){
      hdma.pre[cc] =1
    }
    else{
      hdma.pre[cc] <- tp/(length(hdma.cla))
    }
  }

  #####Multimed
  set.seed(18)
  multm <- try(medTest(E=as.matrix(t),M=as.matrix(log(M_mat+1)),Y=Y_mat,Z=NULL,  w=1),silent = TRUE)
  if(inherits(multm,"try-error")){multm.rec[cc]=multm.pre[cc]=NA}else{
    multm.cla <- which(p.adjust(multm[,2],"fdr") <0.05)
    tp <- length(multm.cla[multm.cla%in%dif])
    fp <- length(multm.cla)-tp
    fn <- length(dif)-tp
    multm.rec[cc] <- tp/length(dif)
    if(length(multm.cla)==0 ){
      multm.pre[cc] =1
    }
    else{
      multm.pre[cc] <- tp/(length(multm.cla))
    }
  }







  #####LDM
  source("./LDM_fun.R")
  meta <-data.frame(trt=trt,Y=Y)
  ldm.res <- ldm(M_mat~ trt+Y,
                 data=meta, seed=67817,fdr.nominal=0.05,
                 test.mediation=TRUE)

  ldm.cla <- ldm.res$med.detected.otu.omni
  tp <- length(ldm.cla[ldm.cla%in%dif])
  fp <- length(ldm.cla)-tp
  fn <- length(dif)-tp
  ldm.rec[cc] <- tp/length(dif)
  if(length(ldm.cla)==0){
    ldm.pre[cc] =1
  }
  else{
    ldm.pre[cc] <- tp/(length(dec.cla))
  }





}

write.table(rbind(cramed.rec,cramed.pre,
                  ikt.rec,ikt.pre,hima.rec,hima.pre,
                  hima2.rec,hima2.pre,hdma.rec,hdma.pre,
                  multm.rec,multm.pre,
                  ldm.rec,ldm.pre), "./figure3_b.txt")


dif <- 1:5
cramed.rec=cramed.pre=c()
ikt.rec=ikt.pre=hima.rec=hima.pre=c()
hima2.rec=hima2.pre=hdma.rec=hdma.pre=c()
multm.rec=multm.pre=ldm.rec=ldm.pre=c()
for(cc in 1:100){
  #######generate simulated data
  generate_hnb <- sim_hnb(otu_n=1000, num=200, alpha=-5, beta=5, gamma=-5, conf.num=3)
  M_mat <- generate_hnb$M_mat
  Y <-generate_hnb$Y
  trt <-generate_hnb$trt
  #######CRAmed
  results.cra <- try(CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05, n.times=100,
                            method="fdr",CI=FALSE),silent = TRUE)

  dec.cra <-  results.cra$Mediators
  tp <- length(dec.cra[dec.cra%in%dif])
  fp <- length(dec.cra)-tp
  fn <- length(dif)-tp
  cramed.rec[cc] <- tp/length(dif)
  if(length(dec.cra)==0){
    cramed.pre[cc] =1
  }
  else{
    cramed.pre[cc] <- tp/(length(dec.cra))
  }
  #####IKT
  med.mat = matrix(NA,otu_n,2)
  for(j in 1:otu_n){
    Y_data <- as.data.frame(cbind(trt,log(M_mat[,j]+1),Y))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    b <-  lm(Mediator ~ Treat, data = Y_data)
    c <- lm(Outcom ~ Treat + Mediator, data=Y_data)
    ikt.res <- try(mediate(b, c, sims=100, treat="Treat", mediator="Mediator"), silent = TRUE)
    if(inherits(contcont, "try-error")) {med.mat[j,]=NA}
    else{
      ikt.res <- summary(ikt.res)
      med.mat[j,] <- c(ikt.res$d0.p,ikt.res$z0.p)
    }
    print(j)
  }

  ikt.cla <- which(p.adjust(med.mat[,1],"fdr")<0.05)
  tp <- length(ikt.cla[ikt.cla%in%dif])
  fp <- length(ikt.cla)-tp
  fn <- length(dif)-tp
  ikt.rec[cc] <- tp/length(dif)
  if(length(ikt.cla)==0 ){
    ikt.pre[cc] =1
  }
  else{
    ikt.pre[cc] <- tp/(length(ikt.cla))
  }

  #####HIMA
  M_mat <- as.data.frame(M_mat)
  colnames(M_mat) = seq(1:otu_n)
  hima.fit <- try(hima(X =trt, Y = Y, M =log(M_mat+1), COV.XM =NULL,scale = FALSE,M.family = c( "gaussian"),Y.family = c("gaussian"),
                       verbose = FALSE),silent=TRUE)
  if(inherits(hima.fit,"try-error")){hima.rec[cc]=hima.pre[cc]=NA}else{
    hima.sel <- hima.fit[which(hima.fit$BH.FDR<0.05),]
    hima.cla <- unlist(lapply(strsplit(rownames(hima.sel),'`'),function(x){x[2]}))

    tp <- length(hima.cla[hima.cla%in%dif])
    fp <- length(hima.cla)-tp
    fn <- length(dif)-tp
    hima.rec[cc] <- tp/length(dif)
    if(length(hima.cla)==0 ){
      hima.pre[cc] =1
    }
    else{
      hima.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HIMA2

  source("./HIMA2.R")
  hima2.fit <- try(HIMA2(X=as.matrix(trt),Y=Y,M=log(M_mat+1),Z=NULL,cut=0.05),silent = TRUE)
  if(inherits(hima2.fit,"try-error")){hima2.rec[cc]=hima2.pre[cc]=NA}
  else{
    hima2.cla <- unlist(lapply(strsplit(hima2.fit$M,'M'),function(x){x[2]}))
    tp <- length(hima2.cla[hima2.cla%in%dif])
    fp <- length(hima2.cla)-tp
    fn <- length(dif)-tp
    hima2.rec[cc] <- tp/length(dif)
    if(length(hima2.cla)==0 ){
      hima2.pre[cc] =1
    }
    else{
      hima2.pre[cc] <- tp/(length(hima.cla))
    }
  }

  #####HDMA
  source("./utils.R")
  source("./hdma.R")

  hdma.res <-  try(hdma(trt, Y=Y, M=log(M_mat+1), COV.XM = NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN =NULL,
                        parallel = FALSE, ncore = 1, verbose = FALSE),silent=TRUE)
  if(inherits(hdma.res,"try-error")){
    hdma.rec[cc]=hdma.pre[cc]=NA
  }
  else{
    hdma.cla <- rownames(hdma.res)[which(p.adjust(hdma.res$P.value,"fdr")<0.05)]

    tp <- length(hdma.cla[hdma.cla%in%dif])
    fp <- length(hdma.cla)-tp
    fn <- length(dif)-tp
    hdma.rec[cc] <- tp/length(dif)
    if(length(hdma.cla)==0 ){
      hdma.pre[cc] =1
    }
    else{
      hdma.pre[cc] <- tp/(length(hdma.cla))
    }
  }

  #####Multimed
  set.seed(18)
  multm <- try(medTest(E=as.matrix(t),M=as.matrix(log(M_mat+1)),Y=Y_mat,Z=NULL,  w=1),silent = TRUE)
  if(inherits(multm,"try-error")){multm.rec[cc]=multm.pre[cc]=NA}else{
    multm.cla <- which(p.adjust(multm[,2],"fdr") <0.05)
    tp <- length(multm.cla[multm.cla%in%dif])
    fp <- length(multm.cla)-tp
    fn <- length(dif)-tp
    multm.rec[cc] <- tp/length(dif)
    if(length(multm.cla)==0 ){
      multm.pre[cc] =1
    }
    else{
      multm.pre[cc] <- tp/(length(multm.cla))
    }
  }







  #####LDM
  source("./LDM_fun.R")
  meta <-data.frame(trt=trt,Y=Y)
  ldm.res <- ldm(M_mat~ trt+Y,
                 data=meta, seed=67817,fdr.nominal=0.05,
                 test.mediation=TRUE)

  ldm.cla <- ldm.res$med.detected.otu.omni
  tp <- length(ldm.cla[ldm.cla%in%dif])
  fp <- length(ldm.cla)-tp
  fn <- length(dif)-tp
  ldm.rec[cc] <- tp/length(dif)
  if(length(ldm.cla)==0){
    ldm.pre[cc] =1
  }
  else{
    ldm.pre[cc] <- tp/(length(dec.cla))
  }





}

write.table(rbind(cramed.rec,cramed.pre,
                  ikt.rec,ikt.pre,hima.rec,hima.pre,
                  hima2.rec,hima2.pre,hdma.rec,hdma.pre,
                  multm.rec,multm.pre,
                  ldm.rec,ldm.pre), "./figure3_c.txt")


######figure
library(ggplot2)
library(Rmisc)
library(tidyquant)
aa= read.table("./figure3_a.txt")

f1 <- 2*aa[1+jk,]*aa[2+jk,]/(aa[1+jk,]+aa[2+jk,])
f2 <- 2*aa[3+jk,]*aa[4+jk,]/(aa[3+jk,]+aa[4+jk,])
f3 <- 2*aa[5+jk,]*aa[6+jk,]/(aa[5+jk,]+aa[6+jk,])
f4 <- 2*aa[7+jk,]*aa[8+jk,]/(aa[7+jk,]+aa[8+jk,])
f5 <- 2*aa[9+jk,]*aa[10+jk,]/(aa[9+jk,]+aa[10+jk,])
f6 <- 2*aa[11+jk,]*aa[12+jk,]/(aa[11+jk,]+aa[12+jk,])
f7 <- 2*aa[13+jk,]*aa[14+jk,]/(aa[13+jk,]+aa[14+jk,])
f8 <- 2*aa[15+jk,]*aa[16+jk,]/(aa[15+jk,]+aa[16+jk,])


data1 <-data.frame(value=c(unlist(aa[1+jk,]),
                           unlist(aa[3+jk,]),
                           unlist(aa[5+jk,]),
                           unlist(aa[7+jk,]),
                           unlist(aa[9+jk,]),
                           unlist(aa[11+jk,]),
                           unlist(aa[13+jk,]),
                           unlist(aa[15+jk,])
),method = rep("Recall",800),
group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2+jk,]),
  unlist(aa[4+jk,]),unlist(aa[6+jk,]),
  unlist(aa[8+jk,]),
  unlist(aa[10+jk,]),
  unlist(aa[12+jk,]),
  unlist(aa[14+jk,]),
  unlist(aa[16+jk,])),method = rep("Precision",800),
  group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2),unlist(f3),unlist(f4),unlist(f5),unlist(f6),unlist(f7),unlist(f8)),
                   method = rep("F1",800),
                   group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))


final.data1 <- rbind(data1,data2,data3)

aa= read.table("./figure3_b.txt")

f1 <- 2*aa[1+jk,]*aa[2+jk,]/(aa[1+jk,]+aa[2+jk,])
f2 <- 2*aa[3+jk,]*aa[4+jk,]/(aa[3+jk,]+aa[4+jk,])
f3 <- 2*aa[5+jk,]*aa[6+jk,]/(aa[5+jk,]+aa[6+jk,])
f4 <- 2*aa[7+jk,]*aa[8+jk,]/(aa[7+jk,]+aa[8+jk,])
f5 <- 2*aa[9+jk,]*aa[10+jk,]/(aa[9+jk,]+aa[10+jk,])
f6 <- 2*aa[11+jk,]*aa[12+jk,]/(aa[11+jk,]+aa[12+jk,])
f7 <- 2*aa[13+jk,]*aa[14+jk,]/(aa[13+jk,]+aa[14+jk,])
f8 <- 2*aa[15+jk,]*aa[16+jk,]/(aa[15+jk,]+aa[16+jk,])


data1 <-data.frame(value=c(unlist(aa[1+jk,]),
                           unlist(aa[3+jk,]),
                           unlist(aa[5+jk,]),
                           unlist(aa[7+jk,]),
                           unlist(aa[9+jk,]),
                           unlist(aa[11+jk,]),
                           unlist(aa[13+jk,]),
                           unlist(aa[15+jk,])
),method = rep("Recall",800),
group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2+jk,]),
  unlist(aa[4+jk,]),unlist(aa[6+jk,]),
  unlist(aa[8+jk,]),
  unlist(aa[10+jk,]),
  unlist(aa[12+jk,]),
  unlist(aa[14+jk,]),
  unlist(aa[16+jk,])),method = rep("Precision",800),
  group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2),unlist(f3),unlist(f4),unlist(f5),unlist(f6),unlist(f7),unlist(f8)),
                   method = rep("F1",800),
                   group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))


final.data2 <- rbind(data1,data2,data3)


aa= read.table("./figure3_c.txt")
f1 <- 2*aa[1+jk,]*aa[2+jk,]/(aa[1+jk,]+aa[2+jk,])
f2 <- 2*aa[3+jk,]*aa[4+jk,]/(aa[3+jk,]+aa[4+jk,])
f3 <- 2*aa[5+jk,]*aa[6+jk,]/(aa[5+jk,]+aa[6+jk,])
f4 <- 2*aa[7+jk,]*aa[8+jk,]/(aa[7+jk,]+aa[8+jk,])
f5 <- 2*aa[9+jk,]*aa[10+jk,]/(aa[9+jk,]+aa[10+jk,])
f6 <- 2*aa[11+jk,]*aa[12+jk,]/(aa[11+jk,]+aa[12+jk,])
f7 <- 2*aa[13+jk,]*aa[14+jk,]/(aa[13+jk,]+aa[14+jk,])
f8 <- 2*aa[15+jk,]*aa[16+jk,]/(aa[15+jk,]+aa[16+jk,])


data1 <-data.frame(value=c(unlist(aa[1+jk,]),
                           unlist(aa[3+jk,]),
                           unlist(aa[5+jk,]),
                           unlist(aa[7+jk,]),
                           unlist(aa[9+jk,]),
                           unlist(aa[11+jk,]),
                           unlist(aa[13+jk,]),
                           unlist(aa[15+jk,])
),method = rep("Recall",800),
group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2+jk,]),
  unlist(aa[4+jk,]),unlist(aa[6+jk,]),
  unlist(aa[8+jk,]),
  unlist(aa[10+jk,]),
  unlist(aa[12+jk,]),
  unlist(aa[14+jk,]),
  unlist(aa[16+jk,])),method = rep("Precision",800),
  group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2),unlist(f3),unlist(f4),unlist(f5),unlist(f6),unlist(f7),unlist(f8)),
                   method = rep("F1",800),
                   group=rep(c("CRAmed","Naive-CRAmed","IKT","HIMA","HIMA2","HDMA","MultiMed","LDM-med"),each=100))


final.data3 <- rbind(data1,data2,data3)



final.data <- rbind(final.data1,final.data2,final.data3)

sam.si <- c(rep("ZINB",nrow(final.data1)),rep("Hurdle-Poisson",nrow(final.data2)),
            rep("Hurdle-NB",nrow(final.data3)))


final.data$sam.si <- factor(sam.si,levels=c("ZINB","Hurdle-Poisson","Hurdle-NB"))


fid1 <- summarySE(final.data, measurevar="value", groupvars=c("method","group","sam.si"),na.rm = TRUE)

fid1$method=factor(fid1$method,levels=c("Recall","Precision","F1"))

Method=factor(fid1$group,levels=c("IKT","MultiMed","HIMA","HIMA2","HDMA","LDM-med","Naive-CRAmed","CRAmed"))




library("tidyquant")
theme_set(theme_bw())
pic=ggplot(fid1,aes(x=method,y=value,fill=Method,col=Method))+
  geom_bar(stat = "identity", position=position_dodge(0.6), width=0.6)+#, position=position_dodge(), width=0.6)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1,position=position_dodge(0.6)) +
  ylab("Value") +
  xlab("")+
  scale_y_continuous(breaks=seq(0,1,0.3))+
scale_fill_manual(values =c("#8DD3C7", "#BC80BD","#80B1D3" ,"#FDB462","#A6CEE3","#BEBADA","#F4CAE4",  "#FB8072"
))+
  scale_color_manual(values =c("#8DD3C7", "#BC80BD","#80B1D3" ,"#FDB462","#A6CEE3","#BEBADA","#F4CAE4",  "#FB8072"
  ))+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))+
  theme(axis.text.x =element_text(size=16), axis.text.y=element_text(size=16))+
  theme(legend.text= element_text(size=11))+theme(plot.title = element_text(hjust = 0.5))+#theme(legend.position = "none")+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(legend.key.size=unit(0.7,'cm'))+
  theme(legend.key.width=unit(0.8,'cm'))+
  facet_grid(sam.si~.)
pic
