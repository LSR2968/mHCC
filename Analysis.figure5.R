library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(dplyr)
library(stringr)
library(textshape)
library(tidyr)
library(crayon)
rmOutlier=function(x){
  tmp.sd=sd(x,na.rm=T)
  tmp.rm=c(x[x<median(x,na.rm = T)-tmp.sd], x[x>median(x,na.rm = T)+tmp.sd])
  x[x %in% tmp.rm]=NA
  return(x)
}

data_meta=read.table("Input/Input.meta.txt",header = T,stringsAsFactors = F,check.names = F)
data_bacteria=read.table("Input/尚儒分析用taxa.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)

# re-arrange ID
data_bacteria=data_bacteria[rownames(data_bacteria) %in% data_meta$`sample.id`,]
data_meta=data_meta[data_meta$sample.id %in% rownames(data_bacteria),]
data_meta=data_meta[order(data_meta$`sample.id`),]
data_bacteria=data_bacteria[order(rownames(data_bacteria)),]
rownames(data_bacteria)=data_meta$sample

# phenotype association
result_individual = foreach(i=1:ncol(data_bacteria),.combine = rbind) %do%  {
  tmp.bacteria=colnames(data_bacteria)[i]
  tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],data_meta[,c("patient","sample","ALT","AST","AFP","VI","Tbilirubin","Cbilirubin","P5266","age","sex","BMI","smoking","drinking")],by.x="row.names",by.y="sample",all=F)
  
  cat(tmp.bacteria,"\n")
  colnames(tmp.data)[2]="Bacteria"
  #tmp.data=tmp.data[tmp.data$Row.names %like% "T",]
  tmp.data = foreach(n=1:length(unique(tmp.data$patient)),.combine = rbind) %do%  {
    tmp.patient=unique(tmp.data$patient)[n]
    tmp.bac=mean(tmp.data$Bacteria[tmp.data$patient==tmp.patient])
    tmp.mean=tmp.data[tmp.data$patient==tmp.patient,]
    tmp.mean=tmp.mean[1,]
    tmp.mean$Bacteria=tmp.bac
    return.string=tmp.mean
  }
  
  tmp.data$Bacteria[tmp.data$Bacteria==0]=NA
  tmp.data$Bacteria=log2(tmp.data$Bacteria) # log transform
  tmp.data$Bacteria=rmOutlier(tmp.data$Bacteria)
  tmp.data$P5266[tmp.data$P5266==9]=NA
  tmp.data$ALT=log2(tmp.data$ALT)
  tmp.data$AST=log2(tmp.data$AST)
  tmp.data$AFP=log2(tmp.data$AFP)
  tmp.data$Tbilirubin=log2(tmp.data$Tbilirubin)
  tmp.data$Cbilirubin=log2(tmp.data$Cbilirubin)
  tmp.data[tmp.data==-Inf]=NA
  
  mod1=lm(Bacteria ~ ALT + age + sex + BMI + smoking + drinking,data = tmp.data)
  mod1=as.data.frame(summary(mod1)$coef)
  
  mod2=lm(Bacteria ~ AST + age + sex + BMI + smoking + drinking,data = tmp.data)
  mod2=as.data.frame(summary(mod2)$coef)
  
  mod3=lm(Bacteria ~ AFP + age + sex + BMI + smoking + drinking,data = tmp.data)
  mod3=as.data.frame(summary(mod3)$coef)
  
  mod4=lm(Bacteria ~ Tbilirubin + age + sex + BMI + smoking + drinking,data = tmp.data)
  mod4=as.data.frame(summary(mod4)$coef)
  
  mod5=lm(Bacteria ~ Cbilirubin + age + sex + BMI + smoking + drinking,data = tmp.data)
  mod5=as.data.frame(summary(mod5)$coef)
  
  #mod6=lm(Bacteria ~ P5266 + age + sex + BMI + smoking + drinking,data = tmp.data)
  #mod6=as.data.frame(summary(mod6)$coef)
  
  return.string=data.frame(Bacteria=tmp.bacteria,beta.ALT=mod1$Estimate[rownames(mod1)=="ALT"],tvalue.ALT=mod1$`t value`[rownames(mod1)=="ALT"],pvalue.ALT=mod1$`Pr(>|t|)`[rownames(mod1)=="ALT"],
                           beta.AST=mod2$Estimate[rownames(mod2)=="AST"],tvalue.AST=mod2$`t value`[rownames(mod2)=="AST"],pvalue.AST=mod2$`Pr(>|t|)`[rownames(mod2)=="AST"],
                           beta.AFP=mod3$Estimate[rownames(mod3)=="AFP"],tvalue.AFP=mod3$`t value`[rownames(mod3)=="AFP"],pvalue.AFP=mod3$`Pr(>|t|)`[rownames(mod3)=="AFP"],
                           beta.Tbilirubin=mod4$Estimate[rownames(mod4)=="Tbilirubin"],tvalue.Tbilirubin=mod4$`t value`[rownames(mod4)=="Tbilirubin"],pvalue.Tbilirubin=mod4$`Pr(>|t|)`[rownames(mod4)=="Tbilirubin"],
                           beta.Cbilirubin=mod5$Estimate[rownames(mod5)=="Cbilirubin"],tvalue.Cbilirubin=mod5$`t value`[rownames(mod5)=="Cbilirubin"],pvalue.Cbilirubin=mod5$`Pr(>|t|)`[rownames(mod5)=="Cbilirubin"]
                           #beta.P5266=mod6$Estimate[rownames(mod6)=="P5266"],pvalue.P5266=mod6$`Pr(>|t|)`[rownames(mod6)=="P5266"]
                           )
}

result_individual=result_individual[result_individual$Bacteria %in% c("Bacteroides","Lawsonella","Corynebacterium","Streptococcus",
                                      "Chryseobacterium","Acinetobacter","Muribaculaceae","Sphingomonas","Pseudomonas","Enterococcus"),]

tmp.alt=result_individual[,c("Bacteria","beta.ALT","tvalue.ALT")]
tmp.ast=result_individual[,c("Bacteria","beta.AST","tvalue.AST")]
tmp.afp=result_individual[,c("Bacteria","beta.AFP","tvalue.AFP")]
colnames(tmp.alt)=c("Bacteria","Beta","tvalue")
colnames(tmp.ast)=c("Bacteria","Beta","tvalue")
colnames(tmp.afp)=c("Bacteria","Beta","tvalue")
tmp.alt$Cate="ALT"
tmp.ast$Cate="AST"
tmp.afp$Cate="AFP"
tmp.plot=rbind(tmp.afp,tmp.alt,tmp.ast)
tmp.plot$Bacteria=factor(tmp.plot$Bacteria,levels = rev(c("Bacteroides","Lawsonella","Corynebacterium","Streptococcus",
                                                      "Chryseobacterium","Acinetobacter","Muribaculaceae","Sphingomonas","Pseudomonas","Enterococcus")))

ggplot(tmp.plot, aes(x = Bacteria, y = tvalue)) +
  geom_segment(
    aes(x = Bacteria, xend = Bacteria, y = 0, yend = tvalue), 
    color = "grey50"
  ) + 
  geom_point(size = 3,aes(color=Cate))  +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )+
  facet_grid(. ~ Cate, scales = "free_y", space = "free_y")+geom_hline(yintercept = 0)+scale_color_jama()+coord_flip()+guides(color=F)
ggsave("OutputFigure/Segment.plot.pdf",width = 10,height = 3)

ggplot(tmp.data, aes(P5266, Bacteria)) +
  geom_jitter(shape = 21, 
              color = "black", size = 2)+
  geom_smooth(method = lm)+scale_color_jama()+
  theme_bw()

tmp.data=tmp.data[!is.na(tmp.data$P5266),]
tmp.data=tmp.data[(tmp.data$P5266!=9),]
tmp.data$P5266=as.factor(tmp.data$P5266)
tmp.data$Bacteria=rmOutlier(tmp.data$Bacteria)
ggplot(tmp.data, aes(x=P5266, y=Bacteria,fill=P5266)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  xlab("")+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')
wilcox.test(tmp.data$Bacteria[tmp.data$P5266=="0"],tmp.data$Bacteria[tmp.data$P5266=="1"])
mod=lm(Bacteria ~ P5266 + age + sex + BMI + smoking + drinking,data = tmp.data)
summary(mod)

































