library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(phyloseq)
library(dplyr)
source("Microbiome.function.R")
library(nlme)
library(lefser)
require(lme4)

rmOutlier=function(x){
  tmp.sd=sd(x,na.rm=T)
  tmp.rm=c(x[x<median(x)-tmp.sd,x>median(x)+tmp.sd])
  x[x %in% tmp.rm]=NA
  return(x)
}

data_bacteria=read.table("Input/尚儒分析用taxa.txt",row.names = 1,sep = "\t",stringsAsFactors = F,check.names = F,header = T)

# re-arrange ID
data_meta=read.table("Input/Input.meta.txt",stringsAsFactors = F,header = T)
data_bacteria=data_bacteria[rownames(data_bacteria) %in% data_meta$`sample.id`,]
data_meta=data_meta[data_meta$sample.id %in% rownames(data_bacteria),]
data_meta=data_meta[order(data_meta$`sample.id`),]
data_bacteria=data_bacteria[order(rownames(data_bacteria)),]
rownames(data_bacteria)=data_meta$sample

# comparison
compare_all_H=foreach(i=1:ncol(data_bacteria),.combine = rbind) %do%  {
  tmp.bacteria=colnames(data_bacteria)[i]
  tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],data_meta[,c("patient","sample","Group")],by.x="row.names",by.y="sample",all=F)
  colnames(tmp.data)[2]="Bacteria"
  tmp.data$Bacteria=log2(tmp.data$Bacteria) # log transformation
  tmp.data$Bacteria[tmp.data$Bacteria==-Inf]=NA
  tmp.data=tmp.data[!is.na(tmp.data$Bacteria),]
  tmp.data$Bacteria=rmOutlier(tmp.data$Bacteria) # optional
  tmp.data=tmp.data[!is.na(tmp.data$Bacteria),]
  
  if(length(unique(tmp.data$Group))==3){
    mod1=(lme(Bacteria ~ Group, random = ~ 1|patient + age +sex + BMI,tmp.data[tmp.data$Group!="N",])) # H is reference
    mod2=(lme(Bacteria ~ Group, random = ~ 1|patient + age +sex + BMI,tmp.data[tmp.data$Group!="T",])) # H is reference
    
    ano1=anova(mod1)
    mod1=summary(mod1)
    ano2=anova(mod2)
    mod2=summary(mod2)
    
    cat(tmp.bacteria,"\n")
    return.string=data.frame(Bacteria=tmp.bacteria,Tumor.beta=mod1$tTable[,1][2],Tumor.se=mod1$tTable[,2][2],Tumor.Pvalue=ano1$`p-value`[2],
                             Peri.beta=mod2$tTable[,1][2],Peri.se=mod2$tTable[,2][2],Peri.Pvalue=ano2$`p-value`[2])
  }else{return.string=data.frame(Bacteria=tmp.bacteria,Tumor.beta=NA,Tumor.se=NA,Tumor.Pvalue=NA,
                                 Peri.beta=NA,Peri.se=NA,Peri.Pvalue=NA)}

}
compare_all_H$FDR.tumor=p.adjust(compare_all_H$Tumor.Pvalue,method = "BH")
compare_all_H$FDR.peri=p.adjust(compare_all_H$Peri.Pvalue,method = "BH")
write.table(compare_all_H,file = "OutputTable /Compare.bacteria.txt",sep = "\t",row.names = F,quote = F)

compare_T=compare_all_H[,c("Bacteria","Tumor.beta","Tumor.se","Tumor.Pvalue")]
compare_N=compare_all_H[,c("Bacteria","Peri.beta","Peri.se","Peri.Pvalue")]
colnames(compare_N)=c("Bacteria","Beta","SE","Pvalue")
colnames(compare_T)=c("Bacteria","Beta","SE","Pvalue")
compare_N$Up=compare_N$Beta+1.96*compare_N$SE
compare_N$Down=compare_N$Beta-1.96*compare_N$SE
compare_T$Up=compare_T$Beta+1.96*compare_T$SE
compare_T$Down=compare_T$Beta-1.96*compare_T$SE
compare_N$Group="N vs. H"
compare_T$Group="T vs. H"
compare_plot=rbind(compare_N,compare_T)
compare_plot=na.omit(compare_plot)

keep=compare_all_H$Bacteria[compare_all_H$Peri.Pvalue<0.05 | compare_all_H$Tumor.Pvalue<0.05]
keep=na.omit(keep)
compare_plot=compare_plot[compare_plot$Bacteria %in% keep,]
compare_plot=compare_plot[compare_plot$Bacteria!="Unnown",]
compare_T=compare_T[compare_T$Bacteria %in% compare_plot$Bacteria,]
compare_T=compare_T[order(compare_T$Beta,decreasing = T),]
compare_plot$Bacteria=factor(compare_plot$Bacteria,levels = compare_T$Bacteria)
ggplot(compare_plot, aes(Bacteria, Beta)) +
  geom_pointrange(
    aes(ymin = Up, ymax = Down, color = Group),
    position = position_dodge(0.6)
  )+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+guides(fill=F)+coord_flip()+geom_hline(yintercept = 0,linetype="dotted", 
                                                                                                           color = "grey", size=1.5)
ggsave("OutputFigure/Compare.pdf",width = 3,height = 5)

# Lactobacillus, Enterococcus, Muribaculaceae, Bacteroides, only significant in tumor
tmp.bacteria="Lactobacillus"
tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],data_meta[,c("patient","sample","Group")],by.x="row.names",by.y="sample",all=F)
colnames(tmp.data)[2]="Bacteria"
tmp.data$Bacteria=log2(tmp.data$Bacteria)
tmp.data$Bacteria[tmp.data$Bacteria==-Inf]=NA
tmp.data=tmp.data[!is.na(tmp.data$Bacteria),]
tmp.data$Bacteria=rmOutlier(tmp.data$Bacteria)
tmp.data=tmp.data[!is.na(tmp.data$Bacteria),]
ggplot(tmp.data, aes(x=Group, y=Bacteria,fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  xlab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+guides(fill=F)+ylab("Lactobacillus")
ggsave("OutputFigure//Lactobacillus.pdf",width = 2.5,height = 2.5)











































