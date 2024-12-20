library(CRAmed)
library(pscl);packageVersion("pscl")# ‘1.5.5’
# Load data
data(Infant)
# remove missing values
Infant_sub <- subset_samples(Infant, !is.na(weight_growth_pace_during_three_years))
# filter out taxa with a prevalence of less than 10\%
prevdf <- apply(X = otu_table(Infant_sub),
                MARGIN = ifelse(taxa_are_rows(Infant_sub), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(Infant_sub),
                     tax_table(Infant_sub))
prevalenceThreshold = floor(0.1* nsamples(Infant_sub))
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
filter.infant <-  prune_taxa(keepTaxa, Infant_sub)
# metadata of Infant data
meta.infant <- sample_data(filter.infant)
# microbiome data
M_mat <- t((otu_table(filter.infant)@.Data))
# outcome
Y <- as.numeric(meta.infant$weight_growth_pace_during_year_one)
#treatment
trt <- ifelse(meta.infant$csection == TRUE, 1, 0)


####CRAmed
results.infant <- try(CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), CI=TRUE))
saveRDS(results.infant,file="./results.infant.rds")


####Naive.CRAmed
results.infant00 <- try(Naive.CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt)))


####HIMA
library(HIMA)
library(foreach)
colnames(M_mat) = seq(1:ncol(M_mat))

hima.fit <- hima(X=trt, Y=Y, M=log(M_mat@.Data+1), COV.XM = NULL, Y.family = c("gaussian"),M.family = c("gaussian"), method="none",verbose = TRUE,scale=FALSE)
hima.sel <- hima.fit[which(hima.fit$BH.FDR<0.05),]
res.hima.log <- rownames(hima.sel)

####HDMA
source("/Users/tiantian/Documents/causal/hdma.R")
source("/Users/tiantian/Documents/causal/utils.R")
hdma.res <-  try(hdma(X=trt, Y=Y, M=log(M_mat@.Data+1), COV.XM =NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN =NULL,
                      parallel = FALSE, ncore = 1, verbose = FALSE),silent=TRUE)

res.hdma.log <- rownames(hdma.res)[which(p.adjust(hdma.res$P.value,"fdr")<0.05)]

# #########HIMA2
source("/Users/tiantian/Documents/causal/HIMA2.R")

hima2.fit <- try(HIMA2(X=matrix(trt),Y=Y,M=log(M_mat@.Data+1),Z= NULL,cut=0.05),silent = TRUE)
res.hima2.log <- unlist(lapply(strsplit(hima2.fit$M,'M'),function(x){x[2]}))

##############MultiMed
library(MultiMed)
set.seed(18)
multm <- medTest(E=as.matrix(trt),M=as.matrix(log(M_mat@.Data+1)),Y=as.matrix(Y),Z=NULL, w=1)
res.multm.log <- which(p.adjust(multm[,2],"fdr") <0.05)

##############LDM
M_mat=as.matrix(M_mat)
source("./LDM_fun.R")
library(permute)
library(vegan)

meta <-data.frame(trt=trt,Y=Y)
res.ldm.med <- ldm(M_mat~ trt+Y,
                   data=meta, seed=67817,fdr.nominal=0.05,
                   test.mediation=TRUE)

res.ldm <- res.ldm.med$med.detected.otu.omni

##############IKT
library(mediation)
otu_n <- ncol(M_mat)
med.mat = matrix(NA,otu_n,2)
for(j in 1:otu_n){
  set.seed(2023)
  Y_data <- as.data.frame(cbind(trt,log(M_mat@.Data[,j]+1),Y))
  colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
  b <-  lm(Mediator ~ Treat, data = Y_data)
  c <- lm(Outcom ~ Treat + Mediator, data=Y_data)
  ikt.res <- try(mediate(b, c, sims=100, treat="Treat", mediator="Mediator"), silent = TRUE)
  if(inherits(ikt.res, "try-error")) {med.mat[j,]=NA}
  else{
    ikt.res <- summary(ikt.res)
    med.mat[j,] <- c(ikt.res$d0.p,ikt.res$z0.p)
  }
  print(j)
}
res.ikt.log <- which(p.adjust(med.mat[,1],"fdr")<0.05)

dataForUpSetPlot <-list("IKT"=as.character(res.ikt.log),"MultiMed"=as.character(res.multm.log),
                        "HIMA"=as.character(res.hima.log),
                        "HIMA2"=as.character(res.hima2.log),"HDMA"=as.character(res.hdma.log),
                        "LDM-med"=as.character(res.ldm),"CRAmed"=as.character(mediation_set))


colo1=c("#8DD3C7",  "#80B1D3", "#A6CEE3","#FB8072","#BEBADA","#FDB462", "#BC80BD")

library(UpSetR)
library(ggplot2)

pic1=upset(fromList(dataForUpSetPlot),
           nsets=length(dataForUpSetPlot),
           keep.order = TRUE,
           point.size = 1.8,
           line.size = 0.7,
           number.angles = 0,
           queries= list(list(query=intersects,params=list("CRAmed"),color="#FB8072",active=T)),
           text.scale = rep(1.8,6),
           matrix.color="black",
           main.bar.color = 'black',
           mainbar.y.label = 'Intersection of mediators',sets.x.label = "Total number of mediators",
           sets.bar.color=colo1)

pic1


media.infant <- readRDS("./results.infant.rds")

nie.names<- c("4449054:g_Bacteroides","177986:f_Lachnospiraceae","191251:g_Parabacteroides","198511:f_Lachnospiraceae",
              "4442459:g_Parabacteroides","157338:o_Clostridiales","4372511:g__Bacteroides","4356080:f_Barnesiellaceae",
              "846127:g_Parabacteroides","4481613:s_aerofaciens")

sort.names<- c("4449054:g_Bacteroides","191251:g_Parabacteroides","4442459:g_Parabacteroides","4481613:s_aerofaciens",
               "198511:f_Lachnospiraceae","157338:o_Clostridiales","4356080:f_Barnesiellaceae","177986:f_Lachnospiraceae",
               "4372511:g__Bacteroides","846127:g_Parabacteroides")

mediation_set <- media.infant$Mediators
nie.total <- c(media.infant$NIE,media.infant$NIEA,media.infant$NIEP)
nie.ci <- cbind(media.infant$NIE.CI,media.infant$NIEA.CI,media.infant$NIEP.CI)

cramed.nie <-data.frame(NIE=nie.total,method=rep(c("NIE","NIEA","NIEP"),each=length(mediation_set)),lower=nie.ci[1,],upper= nie.ci[2,])
final.data <- data.frame(cramed.nie,microb=nie.names)

final.data$microb <- factor(final.data$microb,levels=sort.names)
final.data$method <- factor(final.data$method,levels=c("NIE","NIEA","NIEP"))


pic2 <- ggplot(final.data,aes(x=microb,y=NIE,col=method))+
  geom_point(alpha = 1,size=2.2)+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(.9),width=.1)+
  theme_classic()+
  xlab(NULL)+ylab('Confidence Interval')+
  scale_color_manual(values=c( "#3C5488B2","#00A087B2",
                               "#F39B7FB2","#91D1C2B2"))+
  scale_fill_manual(values=c( "#3C5488B2","#00A087B2",
                              "#F39B7FB2","#91D1C2B2"))+
  facet_grid(.~method, scales = 'free_x')+
  theme(
    legend.position ="none",
    axis.text = element_text(color='black'),
    #axis.title =
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(face="italic",size = 10),axis.title =element_text(size =14))+
  theme(axis.text.x =element_text(size = 12,angle = 90,hjust = 1))+
  theme(strip.text.y = element_text(angle = 0))+
  geom_hline(yintercept = 0,lty=2,lwd=0.5)
pic2
