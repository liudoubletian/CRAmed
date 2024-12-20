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
naive.cramed.rec=naive.cramed.pre=c()

#######generate simulated data
for(cc in 1:100){
generate_zinb <- sim_zinb(otu_n=100, num=100, alpha=-2, beta=2, gamma=-2)
M_mat <- generate_zinb$M_mat
Y <-generate_zinb$Y
trt <-generate_zinb$trt
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

results.naive <- try(Naive.CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05),silent = TRUE)
dec.naive <-  results.naive$Mediators

tp <- length(dec.naive[dec.naive%in%dif])
fp <- length(dec.naive)-tp
fn <- length(dif)-tp
naive.cramed.rec[cc] <- tp/length(dif)
if(length(dec.naivedec.naive)==0){
  naive.cramed.pre[cc] =1
}
else{
  naive.cramed.pre[cc] <- tp/(length(dec.cla))
}

}

write.table(rbind(cramed.rec,cramed.pre,
                  naive.cramed.rec,naive.cramed.pre), "./figure1_a.txt")


for(cc in 1:100){
  generate_zinb <- sim_zinb(otu_n=1000, num=100, alpha=-2, beta=2, gamma=-2)
  M_mat <- generate_zinb$M_mat
  Y <-generate_zinb$Y
  trt <-generate_zinb$trt
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

  results.naive <- try(Naive.CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05),silent = TRUE)
  dec.naive <-  results.naive$Mediators

  tp <- length(dec.naive[dec.naive%in%dif])
  fp <- length(dec.naive)-tp
  fn <- length(dif)-tp
  naive.cramed.rec[cc] <- tp/length(dif)
  if(length(dec.naivedec.naive)==0){
    naive.cramed.pre[cc] =1
  }
  else{
    naive.cramed.pre[cc] <- tp/(length(dec.cla))
  }

}

write.table(rbind(cramed.rec,cramed.pre,
                  naive.cramed.rec,naive.cramed.pre), "./figure1_b.txt")


for(cc in 1:100){
  generate_zinb <- sim_zinb(otu_n=100, num=200, alpha=-2, beta=2, gamma=-2)
  M_mat <- generate_zinb$M_mat
  Y <-generate_zinb$Y
  trt <-generate_zinb$trt
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

  results.naive <- try(Naive.CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05),silent = TRUE)
  dec.naive <-  results.naive$Mediators

  tp <- length(dec.naive[dec.naive%in%dif])
  fp <- length(dec.naive)-tp
  fn <- length(dif)-tp
  naive.cramed.rec[cc] <- tp/length(dif)
  if(length(dec.naivedec.naive)==0){
    naive.cramed.pre[cc] =1
  }
  else{
    naive.cramed.pre[cc] <- tp/(length(dec.cla))
  }

}

write.table(rbind(cramed.rec,cramed.pre,
                  naive.cramed.rec,naive.cramed.pre), "./figure1_c.txt")


for(cc in 1:100){
  generate_zinb <- sim_zinb(otu_n=1000, num=200, alpha=-2, beta=2, gamma=-2)
  M_mat <- generate_zinb$M_mat
  Y <-generate_zinb$Y
  trt <-generate_zinb$trt
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

  results.naive <- try(Naive.CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), FDR = 0.05),silent = TRUE)
  dec.naive <-  results.naive$Mediators

  tp <- length(dec.naive[dec.naive%in%dif])
  fp <- length(dec.naive)-tp
  fn <- length(dif)-tp
  naive.cramed.rec[cc] <- tp/length(dif)
  if(length(dec.naivedec.naive)==0){
    naive.cramed.pre[cc] =1
  }
  else{
    naive.cramed.pre[cc] <- tp/(length(dec.cla))
  }

}

write.table(rbind(cramed.rec,cramed.pre,
                  naive.cramed.rec,naive.cramed.pre), "./figure1_d.txt")

###########figure
library("tidyquant")

aa= read.table("./figure1_a.txt")

f1 <- 2*aa[1,]*aa[2,]/(aa[1,]+aa[2,])
f2 <- 2*aa[3,]*aa[4,]/(aa[3,]+aa[4,])

data1 <-data.frame(value=c(unlist(aa[1,]),

                           unlist(aa[3,])
),method = rep("Recall",200),
group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2,]),

  unlist(aa[4,])),method = rep("Precision",200),
  group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2)),
                   method = rep("F1",200),
                   group=rep(c("CRAmed","Naive-CRAmed"),each=100))

final.data1 <- rbind(data1,data2,data3)


aa= read.table("./figure1_b.txt")


f1 <- 2*aa[1,]*aa[2,]/(aa[1,]+aa[2,])
f2 <- 2*aa[3,]*aa[4,]/(aa[3,]+aa[4,])


data1 <-data.frame(value=c(unlist(aa[1,]),

                           unlist(aa[3,])
),method = rep("Recall",200),
group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2,]),

  unlist(aa[4,])),method = rep("Precision",200),
  group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2)),
                   method = rep("F1",200),
                   group=rep(c("CRAmed","Naive-CRAmed"),each=100))




final.data2 <- rbind(data1,data2,data3)


aa= read.table("./figure1_c.txt")

f1 <- 2*aa[1,]*aa[2,]/(aa[1,]+aa[2,])
f2 <- 2*aa[3,]*aa[4,]/(aa[3,]+aa[4,])


data1 <-data.frame(value=c(unlist(aa[1,]),

                           unlist(aa[3,])
),method = rep("Recall",200),
group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2,]),

  unlist(aa[4,])),method = rep("Precision",200),
  group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2)),
                   method = rep("F1",200),
                   group=rep(c("CRAmed","Naive-CRAmed"),each=100))



final.data3 <- rbind(data1,data2,data3)


aa= read.table("./figure1_d.txt")

f1 <- 2*aa[1,]*aa[2,]/(aa[1,]+aa[2,])
f2 <- 2*aa[3,]*aa[4,]/(aa[3,]+aa[4,])


data1 <-data.frame(value=c(unlist(aa[1,]),

                           unlist(aa[3,])
),method = rep("Recall",200),
group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data2 <-data.frame(value=c(
  unlist(aa[2,]),

  unlist(aa[4,])),method = rep("Precision",200),
  group=rep(c("CRAmed","Naive-CRAmed"),each=100))

data3 <-data.frame(value=c(unlist(f1),unlist(f2)),
                   method = rep("F1",200),
                   group=rep(c("CRAmed","Naive-CRAmed"),each=100))


final.data4 <- rbind(data1,data2,data3)

final.data <- rbind(final.data1,final.data2,final.data3,final.data4)

sam.si <- c(rep("n=100",nrow(final.data1)),rep("n=200",nrow(final.data2)),
            rep("n=100",nrow(final.data3)),rep("n=200",nrow(final.data4)))

taxa.si <- rep(c("m=100",
                 "m=1000"),each=nrow(final.data1)*2)

final.data$sam.si <- factor(sam.si,levels=c("n=100","n=200"))

final.data$taxa.si <- factor(taxa.si,levels=c("m=100",
                                              "m=1000"))
fid1 <- summarySE(final.data, measurevar="value", groupvars=c("method","group","sam.si", "taxa.si"),na.rm = TRUE)

fid1$method=factor(fid1$method,levels=c("Recall","Precision","F1"))
fid1$group <-factor(fid1$group,levels=c("CRAmed","Naive-CRAmed"))
fid1$value <- round(fid1$value,2)

theme_set(theme_bw())
pic=ggplot(fid1,aes(x=method,y=value,group= group))+
  geom_bar(stat = "identity", position=position_dodge(0.6), width=0.6,alpha=1,aes(fill=group))+#, position=position_dodge(), width=0.6)+
  ylab("Value") +
  xlab("")+
  ggtitle("ZINB")+
  scale_y_continuous(breaks=seq(0,2,0.3))+
scale_fill_manual(values =c("#FFBE7A","#8ECFC9"
))+
  scale_color_manual(values =c("#FFBE7A","#8ECFC9"
  ))+
  guides(fill = guide_legend(title = 'Method'))+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16))+
  theme(axis.text.x =element_text(size=16), axis.text.y=element_text(size=16))+
  theme(legend.text= element_text(size=11))+theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(legend.key.size=unit(0.7,'cm'))+
  theme(legend.key.width=unit(0.8,'cm'))+#geom_text(aes(label=value), position=position_stack(vjust=0.5)) +
  facet_grid(sam.si ~  taxa.si)
pic

