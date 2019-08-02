## Go through CML Simulation Results

rm(list=ls())

library(fitdistrplus)
library(ggplot2)
library(DescTools)
library(reshape)
library(pROC)

setwd("C:/Users/Scott/Box Sync/MutBias/Resubmit/Figures/Figure4/Figure4DE/")

## Parse through data

alpha.df = read.csv("AlphaValues121618.csv",stringsAsFactors=F,header=T)
vars = unique(alpha.df$genotype)[2:length(unique(alpha.df$genotype))]

drugs = c("Imatinib","Dasatinib","Nilotinib","Bosutinib")
# drugs = c("Imatinib","Nilotinib")
bias = c(T,F)

summary.df = expand.grid(var=vars,bias=bias,drug=drugs)
summary.df = summary.df[,3:1]

CI = 0.95 # confidence interval
Zalph = qnorm(1-(1-CI)/2)
method="sisonglaz" # wald

files = list.files()
for (i in 1:length(drugs)) {
  drug = drugs[i]
  for (j in 1:length(bias)) {
    
    biasstr = ifelse(bias[j],"WBias","WoBias")
    
    filename = grep(paste("LSCMutBias",drug,biasstr,"011619.csv",sep=""),files,value=T)
    simout.df = read.csv(filename,header=F)
    colnames(simout.df) = c("Cmax","Cmin","log.pop.treat","niter","varID")
    
    simout.df$var[simout.df$varID!=0] = vars[simout.df$varID[simout.df$varID!=0]]
    simout.df$var[simout.df$var==0] = NA
    
    prev = table(simout.df$var)
    prev = prev[vars]
    names(prev) = vars
    prev[is.na(prev)] = 0
    
    prev.est = MultinomCI(prev,conf.level=CI,method=method)
    
    df.row = summary.df$drug==drugs[i]&summary.df$bias==bias[j]
    
    summary.df$sim.pres[df.row] = sum(prev,na.rm=T)/nrow(simout.df)
    summary.df$sim.pest[df.row] = prev/nrow(simout.df)
    summary.df$sim.prev[df.row] = prev
    summary.df$sim.rel.prev[df.row] = prev.est[,"est"]
    summary.df$sim.prev.lci[df.row] = prev.est[,"lwr.ci"]
    summary.df$sim.prev.uci[df.row] = prev.est[,"upr.ci"]
    summary.df$sim.rank[df.row] = rank(-prev)
    
  }
}

summary.df[is.na(summary.df)] = 0

## Imatinib Analysis

imat.df = read.csv("Combined_data_frame_IC_Mutprob_abundance.csv")

clin.prev = imat.df$Abundance
names(clin.prev) = imat.df$Compound
clin.prev = clin.prev[vars]

clin.prev.est = MultinomCI(clin.prev,conf.level=CI,method=method)
clin.rel.prev = clin.prev.est[,"est"]
clin.prev.lci = clin.prev.est[,"lwr.ci"]
clin.prev.uci = clin.prev.est[,"upr.ci"]

clin.rank = rank(-clin.prev)

summary.df$clin.prev[summary.df$drug=="Imatinib"] = clin.prev
summary.df$clin.rel.prev[summary.df$drug=="Imatinib"] = clin.rel.prev
summary.df$clin.prev.lci[summary.df$drug=="Imatinib"] = clin.prev.lci
summary.df$clin.prev.uci[summary.df$drug=="Imatinib"] = clin.prev.uci
summary.df$clin.rank[summary.df$drug=="Imatinib"] = clin.rank

## Nilotinib/Dasatinib Analysis
unique(summary.df$sim.pres[summary.df$drug=="Imatinib"&summary.df$bias==T])
unique(summary.df$sim.pres[summary.df$drug=="Imatinib"&summary.df$bias==F])
unique(summary.df$sim.pres[summary.df$drug=="Nilotinib"&summary.df$bias==T])
unique(summary.df$sim.pres[summary.df$drug=="Nilotinib"&summary.df$bias==F])
unique(summary.df$sim.pres[summary.df$drug=="Dasatinib"&summary.df$bias==T])
unique(summary.df$sim.pres[summary.df$drug=="Dasatinib"&summary.df$bias==F])
unique(summary.df$sim.pres[summary.df$drug=="Bosutinib"&summary.df$bias==T])
unique(summary.df$sim.pres[summary.df$drug=="Bosutinib"&summary.df$bias==F])


nildas.df = read.csv("NilotinibDasatinibPrevalenceData.csv")
nildas.df = nildas.df[match(vars,nildas.df$Compound),]

summary.df$clin.prev[summary.df$drug=="Nilotinib"] = nildas.df$Nilotinib.Abundance
summary.df$clin.rel.prev[summary.df$drug=="Nilotinib"] = nildas.df$Nilotinib.Abundance/sum(nildas.df$Nilotinib.Abundance)

summary.df$clin.rank[summary.df$drug=="Nilotinib"] = rank(-nildas.df$Nilotinib.Abundance)

summary.df$clin.prev[summary.df$drug=="Dasatinib"] = nildas.df$Dasatinib.Abundance
summary.df$clin.rel.prev[summary.df$drug=="Dasatinib"] = nildas.df$Dasatinib.Abundance/sum(nildas.df$Dasatinib.Abundance)

summary.df$clin.rank[summary.df$drug=="Dasatinib"] = rank(-nildas.df$Dasatinib.Abundance)


## Imat/Nil/Dasat Absolute Numbers

n.res.imat = sum(summary.df$clin.prev[summary.df$drug=="Imatinib"&summary.df$bias==T])
n.res.nilot = sum(summary.df$clin.prev[summary.df$drug=="Nilotinib"&summary.df$bias==T])
n.res.dasat = sum(summary.df$clin.prev[summary.df$drug=="Dasatinib"&summary.df$bias==T])

summary.df$sim.pred[summary.df$drug=="Imatinib"] = summary.df$sim.rel.prev[summary.df$drug=="Imatinib"]*n.res.imat
summary.df$sim.pred[summary.df$drug=="Nilotinib"] = summary.df$sim.rel.prev[summary.df$drug=="Nilotinib"]*n.res.nilot
summary.df$sim.pred[summary.df$drug=="Dasatinib"] = summary.df$sim.rel.prev[summary.df$drug=="Dasatinib"]*n.res.dasat

abs.lim = max(summary.df[,c("sim.pred","clin.prev")],na.rm=T)

abs.fit.wbias = lm(sim.pred~clin.prev,data=summary.df[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==T,])
m.abs.wbias = abs.fit.wbias$coefficients[2]
b.abs.wbias = abs.fit.wbias$coefficients[1]

abs.fit.wobias = lm(sim.pred~clin.prev,data=summary.df[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==F,])
m.abs.wobias = abs.fit.wobias$coefficients[2]
b.abs.wobias = abs.fit.wobias$coefficients[1]

ggplot(summary.df[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==T,],aes(x=clin.prev,y=sim.pred))+theme_bw()+
  geom_abline(slope=1,color="gray75",size=2,linetype=2)+
  geom_abline(slope=m.abs.wbias,intercept=b.abs.wbias,color="red",size=2)+
  geom_point(aes(color=drug),size=5,alpha=.8)+
  # geom_segment(aes(xend=clin.rel.prev,yend=sim.prev.lci),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.rel.prev,yend=sim.prev.uci),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.prev.lci,yend=sim.rel.prev),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.prev.uci,yend=sim.rel.prev),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  xlab("Prevalence: Clinically Observed")+ylab("Prevalence: Model Predicted")+
  ggtitle(paste("ABL1 Resistance Prevalence\nWith Mutation Bias"))+  #: Pearson rho = ",round(pears.wbias,2),sep=""))+
  scale_color_manual("Drug",values=c("#ff9626","#44e627","#22a4d6"))+
  theme(plot.title=element_text(hjust=0.5,size=20,face="bold"),
        axis.title=element_text(size=18,face="bold"),
        axis.text=element_text(size=16,face="bold",color="black"),
        legend.title.align = 0.5,
        legend.title=element_text(size=14,face="bold",color="black"),
        legend.text=element_text(size=14,color="black")
  )+
  xlim(0,abs.lim)+ylim(0,abs.lim)
ggsave("AbsWBias.png",height=6,width=6.5)

ggplot(summary.df[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==F,],aes(x=clin.prev,y=sim.pred))+theme_bw()+
  geom_abline(slope=1,color="gray75",size=2,linetype=2)+
  geom_abline(slope=m.abs.wobias,intercept=b.abs.wobias,color="red",size=2)+
  geom_point(aes(color=drug),size=5,alpha=.8)+
  # geom_segment(aes(xend=clin.rel.prev,yend=sim.prev.lci),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.rel.prev,yend=sim.prev.uci),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.prev.lci,yend=sim.rel.prev),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  # geom_segment(aes(xend=clin.prev.uci,yend=sim.rel.prev),arrow=arrow(angle=90,length=unit(0.25,"cm")),size=2)+
  xlab("Prevalence: Clinically Observed")+ylab("Prevalence: Model Predicted")+
  ggtitle(paste("ABL1 Resistance Prevalence\nWithout Mutation Bias"))+  #: Pearson rho = ",round(pears.wbias,2),sep=""))+
  scale_color_manual("Drug",values=c("#ff9626","#44e627","#22a4d6"))+
  theme(plot.title=element_text(hjust=0.5,size=20,face="bold"),
        axis.title=element_text(size=18,face="bold"),
        axis.text=element_text(size=16,face="bold",color="black"),
        legend.title.align = 0.5,
        legend.title=element_text(size=14,face="bold",color="black"),
        legend.text=element_text(size=14,color="black")
  )+
  xlim(0,abs.lim)+ylim(0,abs.lim)
ggsave("AbsWoBias.png",height=6,width=6.5)


abs.pears.wbias = cor(summary.df$clin.prev[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==T],
                      summary.df$sim.pred[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==T],method="pearson")
# 0.86

abs.pears.wobias = cor(summary.df$clin.prev[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==F],
                       summary.df$sim.pred[summary.df$drug%in%c("Imatinib","Nilotinib","Dasatinib")&summary.df$bias==F],method="pearson")
# 0.68

## ROC Analysis

# nilotinib - empirical downsampling
nsims = 1e3
qual.score.wbias = rep(NA,nsims)
qual.score.wobias = rep(NA,nsims)

vars = unique(summary.df$var)

drug = "Nilotinib"
obs.bool = summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]>0
present.logical = obs.bool

qwbias.df = as.data.frame(matrix(nrow=nsims,ncol=length(vars)+1))
colnames(qwbias.df) = c("iter",as.character(vars))
qwbias.df$iter = 1:nsims
qwobias.df = qwbias.df

TPR.wo = matrix(nrow=nsims,ncol=length(vars))
FPR.wo = TPR.wo
TPR.w = TPR.wo
FPR.w = TPR.wo

auc.w = rep(NA,nsims)
auc.wo = rep(NA,nsims)

set.seed(10)
for (i in 1:nsims) {
  
  sim.wbias = sample(vars,size=sum(summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]),replace=T,prob=summary.df$sim.prev[summary.df$drug==drug&summary.df$bias==T])
  present.wbias.logical = table(sim.wbias)>0
  qual.score.wbias[i] = mean(present.wbias.logical==present.logical)
  
  sim.wobias = sample(vars,size=sum(summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]),replace=T,prob=summary.df$sim.prev[summary.df$drug==drug&summary.df$bias==F])
  present.wobias.logical = table(sim.wobias)>0
  qual.score.wobias[i] = mean(present.wobias.logical==present.logical)
  
  # qwbias.df[i,2:ncol(qwbias.df)] = present.wbias.logical==present.logical
  # qwobias.df[i,2:ncol(qwbias.df)] = present.wobias.logical==present.logical
  # 
  # dwnsmpl.wbias.df[i,2:ncol(dwnsmpl.wbias.df)] = present.wbias.logical
  # dwnsmpl.wobias.df[i,2:ncol(dwnsmpl.wbias.df)] = present.wobias.logical
  
  roc.w.i = roc(obs.bool,table(sim.wbias))
  
  labels.w = obs.bool[order(table(sim.wbias),decreasing=T)]
  TPR.w[i,]=cumsum(labels.w)/sum(labels.w)
  FPR.w[i,]=cumsum(!labels.w)/sum(!labels.w)
  
  labels.wo = obs.bool[order(table(sim.wobias),decreasing=T)]
  TPR.wo[i,]=cumsum(labels.wo)/sum(labels.wo)
  FPR.wo[i,]=cumsum(!labels.wo)/sum(!labels.wo)
  
  auc.w[i] = roc(obs.bool,table(sim.wbias))$auc
  auc.wo[i] = roc(obs.bool,table(sim.wobias))$auc
  
}

auc.w.mean = mean(auc.w)
auc.wo.mean = mean(auc.wo)

ROC.df = data.frame(TPR.w = c(0,colMeans(TPR.w)),
                    FPR.w = c(0,colMeans(FPR.w)),
                    TPR.wo = c(0,colMeans(TPR.wo)),
                    FPR.wo = c(0,colMeans(FPR.wo)))
ROC.df = data.frame(TPR.w = c(0,(TPR.w)),
                    FPR.w = c(0,(FPR.w)),
                    TPR.wo = c(0,TPR.wo),
                    FPR.wo = c(0,FPR.wo))
ggplot(ROC.df)+theme_bw()+
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color="gray75")+
  geom_line(aes(x=FPR.wo,y=TPR.wo,color="Without"),size=3)+
  geom_line(aes(x=FPR.w,y=TPR.w,color="With"),size=3)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(1,0.5,0))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(name="Mutation Bias",values= c("#204D74","#9AC0D6"))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curve: Nilotinib Model Predictive Power")+
  theme(
    plot.title=element_text(hjust=0.5,size=22,face="bold"),
    axis.title=element_text(size=20,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=22,hjust=0.5),
    legend.text=element_text(size=20),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill="gray85")
  )
# ggsave("NilotinibROC.pdf",width=6,height=6)


qwbias.df = melt(qwbias.df,id="iter")
qwobias.df = melt(qwobias.df,id="iter")

ggplot(qwbias.df,aes(x=variable,y=iter))+theme_bw()+
  geom_raster(aes(fill=value))+
  ggtitle("Nilotinib - With Bias")+
  xlab("Mutant")+ylab("Simulations")+
  scale_fill_manual(name="Predicted",values=c("#f74327","#20c4db"))+
  scale_y_continuous(labels=c())+
  theme(
    plot.title=element_text(hjust=0.5,size=18,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=14,hjust=0.5),
    legend.text=element_text(size=10)
    )

ggplot(qwobias.df,aes(x=variable,y=iter))+theme_bw()+
  geom_raster(aes(fill=value))+
  ggtitle("Nilotinib - Without Bias")+
  xlab("Mutant")+ylab("Simulations")+
  scale_fill_manual(name="Predicted",values=c("#f74327","#20c4db"))+
  scale_y_continuous(labels=c())+
  theme(
    plot.title=element_text(hjust=0.5,size=18,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=14,hjust=0.5),
    legend.text=element_text(size=10)
  )


# dasatinib - empirical downsampling
nsims = 1e3
qual.score.wbias = rep(NA,nsims)
qual.score.wobias = rep(NA,nsims)
rocs = rep(NA,nsims)

vars = unique(summary.df$var)

drug = "Dasatinib"
obs.bool = summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]>0
present.logical = obs.bool

qwbias.df = as.data.frame(matrix(nrow=nsims,ncol=length(vars)+1))
colnames(qwbias.df) = c("iter",as.character(vars))
qwbias.df$iter = 1:nsims
qwobias.df = qwbias.df

dwnsmpl.wbias.df = qwbias.df
dwnsmpl.wobias.df = qwbias.df

# roc.w.df = expand.grid(var=1:length(vars),iter=1:nsims)
# roc.w.df = roc.w.df[,c(2,1)]
# roc.w.df = cbind(roc.w.df,
#                  as.data.frame(matrix(nrow=nrow(roc.w.df),ncol=2)))
# colnames(roc.w.df) = c("iter","var","sens","comp.spec")
# roc.wo.df = roc.w.df

auc.w = rep(NA,nsims)
auc.wo = rep(NA,nsims)

TPR.wo = matrix(nrow=nsims,ncol=length(vars))
FPR.wo = TPR.wo
TPR.w = TPR.wo
FPR.w = TPR.wo

set.seed(10)
for (i in 1:nsims) {
  
  sim.wbias = sample(vars,size=sum(summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]),replace=T,prob=summary.df$sim.prev[summary.df$drug==drug&summary.df$bias==T])
  present.wbias.logical = table(sim.wbias)>0
  qual.score.wbias[i] = mean(present.wbias.logical==present.logical)
  
  sim.wobias = sample(vars,size=sum(summary.df$clin.prev[summary.df$drug==drug&summary.df$bias==T]),replace=T,prob=summary.df$sim.prev[summary.df$drug==drug&summary.df$bias==F])
  present.wobias.logical = table(sim.wobias)>0
  qual.score.wobias[i] = mean(present.wobias.logical==present.logical)
  
  qwbias.df[i,2:ncol(qwbias.df)] = present.wbias.logical==present.logical
  qwobias.df[i,2:ncol(qwbias.df)] = present.wobias.logical==present.logical
  
  dwnsmpl.wbias.df[i,2:ncol(dwnsmpl.wbias.df)] = present.wbias.logical
  dwnsmpl.wobias.df[i,2:ncol(dwnsmpl.wbias.df)] = present.wobias.logical
  
  auc.w[i] = roc(obs.bool,table(sim.wbias))$auc
  auc.wo[i] = roc(obs.bool,table(sim.wobias))$auc
  
  labels.w = obs.bool[order(table(sim.wbias),decreasing=T)]
  TPR.w[i,]=cumsum(labels.w)/sum(labels.w)
  FPR.w[i,]=cumsum(!labels.w)/sum(!labels.w)
  
  labels.wo = obs.bool[order(table(sim.wobias),decreasing=T)]
  TPR.wo[i,]=cumsum(labels.wo)/sum(labels.wo)
  FPR.wo[i,]=cumsum(!labels.wo)/sum(!labels.wo)
  
}

auc.w.mean = mean(auc.w)
auc.wo.mean = mean(auc.wo)

ROC.df = data.frame(TPR.w = c(0,colMeans(TPR.w)),
                    FPR.w = c(0,colMeans(FPR.w)),
                    TPR.wo = c(0,colMeans(TPR.wo)),
                    FPR.wo = c(0,colMeans(FPR.wo)))
ggplot(ROC.df)+theme_bw()+
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color="gray75")+
  geom_line(aes(x=FPR.wo,y=TPR.wo,color="Without"),size=3)+
  geom_line(aes(x=FPR.w,y=TPR.w,color="With"),size=3)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c(1,0.5,0))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(name="Mutation Bias",values=c("#307220","#ACD69A"))+
  xlab("Specificity")+ylab("Sensitivity")+
  ggtitle("ROC Curve: Dasatinib Model Predictive Power")+
  theme(
    plot.title=element_text(hjust=0.5,size=22,face="bold"),
    axis.title=element_text(size=20,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=22,hjust=0.5),
    legend.text=element_text(size=20),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill="gray85")
  )
# ggsave("DasatinibROC.pdf",width=6,height=6)


qwbias.df = melt(qwbias.df,id="iter")
qwobias.df = melt(qwobias.df,id="iter")

ggplot(qwbias.df,aes(x=variable,y=iter))+theme_bw()+
  geom_raster(aes(fill=value))+
  ggtitle("Dasatinib - With Bias")+
  xlab("Mutant")+ylab("Simulations")+
  scale_fill_manual(name="Predicted",values=c("#f74327","#20c4db"))+
  scale_y_continuous(labels=c())+
  theme(
    plot.title=element_text(hjust=0.5,size=18,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=14,hjust=0.5),
    legend.text=element_text(size=10)
  )

ggplot(qwobias.df,aes(x=variable,y=iter))+theme_bw()+
  geom_raster(aes(fill=value))+
  ggtitle("Dasatinib - Without Bias")+
  xlab("Mutant")+ylab("Simulations")+
  scale_fill_manual(name="Predicted",values=c("#f74327","#20c4db"))+
  scale_y_continuous(labels=c())+
  theme(
    plot.title=element_text(hjust=0.5,size=18,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    axis.text=element_text(face="bold"),
    legend.title=element_text(size=14,hjust=0.5),
    legend.text=element_text(size=10)
  )


