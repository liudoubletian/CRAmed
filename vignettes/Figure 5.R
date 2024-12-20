library(CRAmed)
library(pscl);packageVersion("pscl")# ‘1.5.5’
# Load data
data(GGMP)
# filter out taxa with a prevalence of less than 10\%
prevdf <- apply(X = otu_table(GGMP),
                MARGIN = ifelse(taxa_are_rows(GGMP), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(GGMP),
                     tax_table(GGMP))
prevalenceThreshold <- floor(0.1* nsamples(GGMP))
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
filter.ggmp <- prune_taxa(keepTaxa, GGMP)
# metadata of GGMP data
meta.ggmp <- sample_data(filter.ggmp)
# microbiome data
M_mat <- t((otu_table(filter.ggmp)@.Data))
# outcome
Y <- as.numeric(unlist(meta.ggmp[,"BMI"]))
# treatment
trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)


####CRAmed
y.names <- c('BMI','Waist','HDL','LDL','TCHO','DBP','TG','FBG','SBP')
for(jj in 1:length(y.names)){
  nam <-y.names[jj]
  Y <- as.numeric(unlist(meta.ggmp[,nam]))
  trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)
  na.y <- sum(is.na(Y))
  M_mat <- t((otu_table(filter.ggmp)@.Data))
  if(na.y==0){
    M_mat <- M_mat
    trt <- trt
    Y <- Y
  }else{
    M_mat <- M_mat[!is.na(Y),]
    trt <- trt[!is.na(Y)]
    Y <- Y[!is.na(Y)]
  }
  set.seed(1)
  samps <- c(which(trt==0),sample(which(trt==1),length(which(trt==0))))
  M_mat <- M_mat[samps,]
  Y <- Y[samps]
  trt <- trt[samps]
  results.ggmp <- try(CRAmed(M_mat=as.matrix(M_mat), Y=Y, Exposure=as.matrix(trt), method="none", CI=FALSE))
  new.name <- paste0("results.ggmp",nam)
  assign(new.name,results.ggmp)
}

ggmp.mediation <- unique(c(results.ggmpBMI$Mediators,results.ggmpWaist$Mediators,results.ggmpDBP$Mediators,results.ggmpTG$Mediators))
#ggmp.mediation <-c(308,581,102,592,801,872)
mediation.names <-tax_table(filter.ggmp)[ggmp.mediation]

p.value=c(0.0002,0.023,0.015,0.003,0.005,0.004,0.010,0.050,0.030,0.040,0.020,0.020,0.020)


library(circlize)
cir.dat <- data.frame('from'=c(rep('Antibiotics',6),'4401580:g_Bacteroides','4436046:g_Dorea',
                               '4425669:g_Coprococcus','4470603:g_Dorea','174749:f_Ruminococcaceae',
                               '4401580:g_Bacteroides','4347159:s_adolescentis'),
                      'to'=c('4436046:g_Dorea',
                             '4425669:g_Coprococcus','4470603:g_Dorea','174749:f_Ruminococcaceae',
                             '4401580:g_Bacteroides','4347159:s_adolescentis',
                             'WC',rep('DBP',3),rep('BMI',2),'TG'),'value'=-log(p.value))

circos.par(gap.after = c("Antibiotics" = 10, "4401580:g_Bacteroides" = 3, "4436046:g_Dorea" = 3, "4425669:g_Coprococcus" = 3, "4470603:g_Dorea" = 3,
                         "174749:f_Ruminococcaceae" = 3, "4401580:g_Bacteroides" = 3, "4347159:s_adolescentis" = 10,
                         "WC" = 3,"DBP" = 3,"BMI" = 3,"TG" = 10))

grid.col <-c('Antibiotics'="#fb8072",'4401580:g_Bacteroides'="#80b1d3",'4436046:g_Dorea'="#80b1d3",
             '4425669:g_Coprococcus'="#80b1d3",'4470603:g_Dorea'="#80b1d3",'174749:f_Ruminococcaceae'="#80b1d3",
             '4347159:s_adolescentis'="#80b1d3",'BMI'="#fdb462",'DBP'="#fdb462",'TG'="#fdb462",
             'WC'="#fdb462")

par(cex = 1.2, mar = c(0, 0, 0, 0))
chordDiagram(cir.dat, directional = 1, direction.type = c( "arrows"),link.arr.length = 0.1,
             grid.col=grid.col)
circos.clear()


# devtools::install_github("ChiLiubio/microeco")
library(microeco)
library(ggplot2)
tax.table <-as.data.frame(tax_table(filter.ggmp))
group <- ifelse(trt==1,"non-antibiotics","antibiotics")
group <- factor(group,levels=c("antibiotics","non-antibiotics"))
sample.table <- data.frame("group"=group,"none"=Y)
rownames(sample.table) <-rownames(M_mat)
M_mat.lefse <-as.data.frame(M_mat[,ggmp.mediation])
lefse.names <-colnames(M_mat)[ggmp.mediation]
tax.table.lefse <- tax.table[lefse.names,]
remove.samps <- which(rowSums(M_mat.lefse)==0)
M_mat.lefse <-as.data.frame(M_mat.lefse[-remove.samps,])
rownames(M_mat.lefse) <-rownames(M_mat)[-remove.samps]
sample.table <- sample.table[-remove.samps,]
colnames(M_mat.lefse) <- c('174749:f_Ruminococcaceae','4401580:g_Bacteroides',
                           '4436046:g_Dorea','4425669:g_Coprococcus',
                           '4470603:g_Dorea','4347159:s_adolescentis')
rownames(tax.table.lefse) <- c('174749:f_Ruminococcaceae','4401580:g_Bacteroides',
                               '4436046:g_Dorea','4425669:g_Coprococcus',
                               '4470603:g_Dorea','4347159:s_adolescentis')
dataset <- microtable$new(sample_table = sample.table,
                          otu_table = as.data.frame(t(M_mat.lefse)),tax_table = tax.table.lefse)
set.seed(1)
lefse <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "group",alpha=1,taxa_level = "b")

color_map <- c("non-antibiotics"="#80B1D3","antibiotics"="#FB8072")
lefseplot <- lefse$plot_diff_bar(use_number = 1:80,
                                 width = 0.8) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()

lefseplot <- lefseplot+scale_color_manual(values=color_map)+scale_fill_manual(values=color_map)+
  theme(legend.position = "bottom",axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),legend.text = element_text(size=15))

lefseplot

library(reshape2)
library(RColorBrewer)
library(ggpubr)

mediation_set1 <-c(308, 581)#BMI
mediation_set2<-c(581)#WC
mediation_set3<-c(102,592,801)#DBP
mediation_set4<-c(872)#TG
nam <-y.names[1]
Y <- as.numeric(unlist(meta.ggmp[,nam]))
trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)
na.y <- sum(is.na(Y))
M_mat <- t((otu_table(filter.ggmp)@.Data))
if(na.y==0){
  M_mat <- M_mat
  trt <- trt
  Y <- Y
}else{
  M_mat <- M_mat[!is.na(Y),]
  trt <- trt[!is.na(Y)]
  Y <- Y[!is.na(Y)]
}
set.seed(1)
samps <- c(which(trt==0),sample(which(trt==1),length(which(trt==0))))
M_mat <- M_mat[samps,]
Y <- Y[samps]
trt <- trt[samps]
line.data <-data.frame(M_mat[,mediation_set1],ymat=Y)
tax_table(filter.ggmp)[mediation_set1]
colnames(line.data) <- c("174749:f_Ruminococcaceae" ,"4401580:g_Bacteroides",'ymat')
line.melt <- melt(line.data, id="ymat")
line.melt$variable <-factor(line.melt$variable,levels=c("4401580:g_Bacteroides","174749:f_Ruminococcaceae"))
colos <- brewer.pal(name="Set1",10)
picline1 <- ggplot(data=line.melt, aes(x=value, y=ymat, group=variable, color=variable)) +
  labs(x=" Observed count", y=nam) +
  geom_point(shape = 1) +
  scale_fill_manual(values =colos )+
  scale_color_manual(values = colos)+
  geom_smooth(method="lm", se = FALSE) +
  theme(panel.grid=element_blank(),  legend.position = "bottom",
        panel.background=element_rect(fill='transparent', color='black'),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
        legend.title = element_text(size=15),legend.text = element_text(size=15)) +
  guides(color=guide_legend(title="Taxa"))
picline1



nam <-y.names[2]
Y <- as.numeric(unlist(meta.ggmp[,nam]))
trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)
na.y <- sum(is.na(Y))
M_mat <- t((otu_table(filter.ggmp)@.Data))
if(na.y==0){
  M_mat <- M_mat
  trt <- trt
  Y <- Y
}else{
  M_mat <- M_mat[!is.na(Y),]
  trt <- trt[!is.na(Y)]
  Y <- Y[!is.na(Y)]
}
set.seed(1)
samps <- c(which(trt==0),sample(which(trt==1),length(which(trt==0))))
M_mat <- M_mat[samps,]
Y <- Y[samps]
trt <- trt[samps]
line.data <-data.frame(M_mat[,mediation_set2],ymat=Y)
tax_table(filter.ggmp)[mediation_set2]
colnames(line.data) <- c("4401580:g_Bacteroides    ",'ymat')
line.melt <- melt(line.data, id="ymat")

colos <- brewer.pal(name="Set1",10)
picline2 <- ggplot(data=line.melt, aes(x=value, y=ymat, group=variable, color=variable)) +
  labs(x=" Observed count", y="WC") +
  geom_point(shape = 1) +
  scale_fill_manual(values =colos )+
  scale_color_manual(values = colos)+
  geom_smooth(method="lm", se = FALSE) +
  theme(panel.grid=element_blank(),  legend.position = "bottom",
        panel.background=element_rect(fill='transparent', color='black'),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
        legend.title = element_text(size=15),legend.text = element_text(size=15)) +
  guides(color=guide_legend(title="Taxa"))
picline2



nam <-y.names[6]
Y <- as.numeric(unlist(meta.ggmp[,nam]))
trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)
na.y <- sum(is.na(Y))
M_mat <- t((otu_table(filter.ggmp)@.Data))
if(na.y==0){
  M_mat <- M_mat
  trt <- trt
  Y <- Y
}else{
  M_mat <- M_mat[!is.na(Y),]
  trt <- trt[!is.na(Y)]
  Y <- Y[!is.na(Y)]
}
set.seed(1)
samps <- c(which(trt==0),sample(which(trt==1),length(which(trt==0))))
M_mat <- M_mat[samps,]
Y <- Y[samps]
trt <- trt[samps]
line.data <-data.frame(M_mat[,mediation_set3],ymat=Y)
tax_table(filter.ggmp)[mediation_set3]
colnames(line.data) <- c( "4436046:g_Dorea","4425669:g_Coprococcus ",  "4470603:g_Dorea" ,'ymat')
line.melt <- melt(line.data, id="ymat")

colos <- c("#E41A1C", "#377EB8", "#984EA3" )
picline3 <- ggplot(data=line.melt, aes(x=value, y=ymat, group=variable, color=variable)) +
  labs(x=" Observed count", y=nam) +
  geom_point(shape = 1) +
  scale_fill_manual(values =colos )+
  scale_color_manual(values = colos)+
  geom_smooth(method="lm", se = FALSE) +
  theme(panel.grid=element_blank(), legend.position = "bottom",
        panel.background=element_rect(fill='transparent', color='black'),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
        legend.title = element_text(size=15),legend.text = element_text(size=15)) +
  guides(color=guide_legend(title="Taxa"))
picline3


nam <-y.names[7]
Y <- as.numeric(unlist(meta.ggmp[,nam]))
trt <- ifelse(meta.ggmp$Antibiotics == "n", 1, 0)
na.y <- sum(is.na(Y))
M_mat <- t((otu_table(filter.ggmp)@.Data))
if(na.y==0){
  M_mat <- M_mat
  trt <- trt
  Y <- Y
}else{
  M_mat <- M_mat[!is.na(Y),]
  trt <- trt[!is.na(Y)]
  Y <- Y[!is.na(Y)]
}
set.seed(1)
samps <- c(which(trt==0),sample(which(trt==1),length(which(trt==0))))
M_mat <- M_mat[samps,]
Y <- Y[samps]
trt <- trt[samps]
line.data <-data.frame(M_mat[,mediation_set4],ymat=Y)
tax_table(filter.ggmp)[mediation_set4]
colnames(line.data) <- c( "4347159:s_adolescentis     "  ,'ymat')
line.melt <- melt(line.data, id="ymat")
colos <- brewer.pal(name="Set1",10)
picline4 <- ggplot(data=line.melt, aes(x=value, y=ymat, group=variable, color=variable)) +
  labs(x=" Observed count", y=nam) +
  geom_point(shape = 1) +
  scale_fill_manual(values =colos )+
  scale_color_manual(values = colos)+
  geom_smooth(method="lm", se = FALSE) +
  theme(panel.grid=element_blank(),  legend.position = "bottom",
        panel.background=element_rect(fill='transparent', color='black'),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),
        legend.title = element_text(size=15),legend.text = element_text(size=15)) +
  guides(color=guide_legend(title="Taxa"))
picline4
ggarrange(picline1,picline2,picline3,picline4,nrow=2,ncol=2)
