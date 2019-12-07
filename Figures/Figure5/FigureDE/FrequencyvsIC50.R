
library(ggplot2)

rm(list=ls())
setwd("C:/Users/Scott/Box Sync/MutBias/Reresubmit/Figures/Figure6/FigureDE/")

gene.vec = c("AR","ESR1","KIT","ALK","EGFR")
col.vec = c("#EC2027","#EC2027","#EC2027","#2E368F","#2E368F")

for (i in 1:length(gene.vec)) {
  
  gene = gene.vec[i]
  
  IC50.df = read.csv(paste(gene,"_IC50Summary.csv",sep=""),header=T,stringsAsFactors=F)
  colnames(IC50.df)[ncol(IC50.df)] = "IC50"
  
  IC50.df = IC50.df[!is.na(IC50.df$IC50),]
  IC50.df$rel.freq = IC50.df$count/sum(IC50.df$count)
  IC50.df$most.common = F
  IC50.df$most.common[which.max(IC50.df$rel.freq)] = T
  
  r = cor(IC50.df$rel.freq,IC50.df$IC50)
  
  ggplot(IC50.df,aes(x=IC50,y=rel.freq))+theme_bw()+
    geom_point(aes(color=most.common),size=4,alpha=0.75)+
    xlab("IC50")+ylab("Frequency")+
    ggtitle(paste(gene,"\nr =",round(r,2)))+
    scale_color_manual(values=c("gray50",col.vec[i]))+
    theme(
      plot.title=element_text(size=20,face="bold",hjust=0.5,color=col.vec[i]),
      axis.title=element_text(size=18,face="bold"),
      axis.text=element_text(size=16,face="bold")
    )+
    guides(color=F)
  
  ggsave(paste(gene,"IC50s.png",sep=""),width=4,height=3)
  
}
