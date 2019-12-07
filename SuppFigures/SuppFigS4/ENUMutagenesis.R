
rm(list=ls())
library(ggplot2)

ENU.df = read.csv("Bradeen.csv",header=T)

# ENU.df = ENU.df[ENU.df$variant%in%our.vars,] # filter for our variants
# rownames(ENU.df) = 1:nrow(ENU.df)

ENU.counts = rep(0,length(our.vars))
names(ENU.counts) = our.vars
for (i in 1:length(ENU.counts)) {
  
  ENU.counts[i] = sum(ENU.df$overall[as.character(ENU.df$variant)==our.vars[i]])
  
}

ENU.plt.df = data.frame(var = our.vars,ENU = ENU.counts,clin = imat.df$Abundance)

ggplot(ENU.plt.df,aes(x=ENU/sum(ENU),y=clin/sum(clin)))+theme_bw()+
  geom_abline(slope=1,color="gray75",size=2,linetype=2)+
  geom_point(aes(color=clin),size=7.5,alpha=0.7,shape=16)+
  xlab("Mutant Frequency:\nMutagenesis Screen")+ylab("Mutant Frequency:\nClinical Population")+
  scale_color_gradient2(low="white",mid="#ff2e00",high="#961c02",midpoint=median(imat.df$Abundance/sum(imat.df$Abundance)))+
  scale_x_continuous(breaks=c(0,0.5),limits=c(0,0.5))+
  scale_y_continuous(breaks=c(0,0.5),limits=c(0,0.5))+
  theme(
    axis.title = element_text(size=28,face="bold"),
    axis.text = element_text(size=22,color="black",face="bold"),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  # xlim(0,0.51)+ylim(0,0.51)+
  guides(color=F)
# ggsave("ENUQuantitative.png",height=7,width=7)

ggplot(ENU.plt.df,aes(x=ENU/sum(ENU),y=clin/sum(clin)))+theme_bw()+
  geom_abline(slope=1,size=2,color="red")+
  geom_point(color='black',size=7.5,alpha=0.7,shape=16)+
  xlab("Mutant Frequency:\nMutagenesis Screen")+ylab("Mutant Frequency:\nClinical Population")+
  scale_color_gradient2(low="white",mid="#ff2e00",high="#961c02",midpoint=median(imat.df$Abundance/sum(imat.df$Abundance)))+
  scale_x_continuous(breaks=c(0,0.5),limits=c(0,0.5))+
  scale_y_continuous(breaks=c(0,0.5),limits=c(0,0.5))+
  theme(
    axis.title = element_text(size=28,face="bold"),
    axis.text = element_text(size=22,color="black",face="bold"),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  # xlim(0,0.51)+ylim(0,0.51)+
  guides(color=F)
