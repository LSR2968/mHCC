
library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(phyloseq)
library(dplyr)
source("Microbiome.function.R")
library(microbiome)
library(crayon)

data_meta=read.table("Input/Input.meta.txt",header = T,stringsAsFactors = F,check.names = F)
data_bacteria=read.table("Input/尚儒分析用taxa.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)

# clr transformation
data_bacteria=as.data.frame(t(data_bacteria))
str(data_bacteria)
data_bacteria=data_bacteria[,colSums(data_bacteria!=0) > (nrow(data_bacteria) * 0.2) ]
data_bacteria=data_bacteria[,!colnames(data_bacteria) %in% c("Other","uncultured")]
data_bacteria=apply(data_bacteria,1,function(x){
  x=x/sum(x)
  return(x)
})
data_bacteria=as.data.frame(t(data_bacteria))
data_bacteria_clr <- zCompositions::cmultRepl(data_bacteria, method="CZM", label=0)
data_bacteria_clr = compositions::clr(data_bacteria_clr)
data_bacteria_clr=as.data.frame(data_bacteria_clr)

# re-arrange ID
data_bacteria=data_bacteria[rownames(data_bacteria) %in% data_meta$`sample.id`,]
data_meta=data_meta[data_meta$sample.id %in% rownames(data_bacteria),]
data_meta=data_meta[order(data_meta$`sample.id`),]
data_bacteria=data_bacteria[order(rownames(data_bacteria)),]
rownames(data_bacteria)=data_meta$sample

# shannon diversity 
detach("package:microbiome", unload=TRUE)
shannon <- as.data.frame(diversity((data_bacteria), index="shannon"))
shannon=merge(shannon,data_meta,by.x="row.names",by.y="sample",all=F)
colnames(shannon)[2]="shannon"
ggplot(shannon, aes(x=Group, y=shannon, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  xlab("")+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') 
ggsave("OutputFigure//Shannono.pdf",width = 3,height = 2.5)
wilcox.test(shannon$shannon[shannon$Group=="T"],shannon$shannon[shannon$Group=="N"])

# beta diversity
beta_diversity=vegdist((data_bacteria),method = "bray")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,data_meta,all = F,by.x="row.names",by.y="sample")
outlier=pca_analysis$Row.names[pca_analysis$V2>5.1 & pca_analysis$Group=="T"]
pca_analysis=pca_analysis[!pca_analysis$Row.names %in% c(outlier,"591L","S1024T2","S945T1","S945T4","S1558T2"),]
ggplot (pca_analysis, aes(V1,V2,color=Group)) + 
  geom_point(size=2) + theme_classic() +
  guides(size=F)+
  scale_color_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  stat_ellipse(type = "norm")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') +xlab("PCA1")+ylab("PCA2")
ggsave("OutputFigure//Bray.pdf",width = 4,height = 3)


# 多灶内部差异
distmat <- as.matrix(beta_diversity)

# cross-patient peri-tumor
tmp.meta=data_meta[data_meta$sample %like% "L",]
dis_cross_peri=distmat[rownames(distmat) %in% tmp.meta$sample,colnames(distmat) %in% tmp.meta$sample]
dis_cross_peri=data.frame(as.table(dis_cross_peri))[lower.tri(dis_cross_peri, diag = F), ]

# cross-patient tumor
tmp.tumor=data_meta[data_meta$sample %like% "T",]
dis_cross_tumor = foreach(i=1:nrow(tmp.tumor),.combine = rbind) %do%  {
  tmp.sample=(tmp.tumor$sample[i])
  tmp.peteint=tmp.tumor$patient[i]
  
  tmp.data=tmp.tumor[!tmp.tumor$patient %in%tmp.peteint,]
  tmp.cross=colnames(distmat)[colnames(distmat) %in% tmp.data$sample]
  tmp.data=distmat[rownames(distmat)==tmp.sample,tmp.cross]
  tmp.data=data.frame(sampleA=tmp.sample,sampleB=tmp.cross,dist=tmp.data)
  rownames(tmp.data)=NULL
  
  return.string=tmp.data
}

# within-patient tumor
tmp.tumor=data_meta[data_meta$sample %like% "T",]
dis_within_tumor = foreach(i=1:nrow(tmp.tumor),.combine = rbind) %do%  {
  tmp.sample=(tmp.tumor$sample[i])
  tmp.peteint=tmp.tumor$patient[i]
  
  tmp.data=tmp.tumor[tmp.tumor$patient %in%tmp.peteint,]
  tmp.cross=colnames(distmat)[colnames(distmat) %in% tmp.data$sample]
  tmp.data=distmat[rownames(distmat)==tmp.sample,tmp.cross]
  tmp.data=data.frame(sampleA=tmp.sample,sampleB=tmp.cross,dist=tmp.data)
  rownames(tmp.data)=NULL
  tmp.data=tmp.data[tmp.data$dist!=0,]
  
  return.string=tmp.data
}

dis_cross_peri$group="Cross-patient N"
dis_cross_tumor$group="Cross-patient T"
dis_within_tumor$group="Within-patient T"
colnames(dis_cross_peri)=c("sampleA","sampleB","dist","group")
distmat=rbind(dis_cross_peri,dis_cross_tumor,dis_within_tumor)

distmat$group=factor(distmat$group,levels = c("Cross-patient N","Cross-patient T","Within-patient T"))
ggplot(distmat, aes(x=group, y=dist,fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  xlab("")+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')
ggsave("OutputFigure//dissimilarity.genus.bray.curtis.pdf",width = 4,height = 3)

wilcox.test(distmat$dist[distmat$group=="Cross-patient peri-tumor"],distmat$dist[distmat$group=="Cross-patient tumor"])
wilcox.test(distmat$dist[distmat$group=="Cross-patient peri-tumor"],distmat$dist[distmat$group=="Within-patient tumor"])
wilcox.test(distmat$dist[distmat$group=="Cross-patient tumor"],distmat$dist[distmat$group=="Within-patient tumor"])

# adonis
ad_taxa=data_bacteria
ad_meta=data_meta[,c("sample","patient","Group","age","sex","BMI","smoking","drinking","Tbilirubin","Cbilirubin","ALT","AST","AFP")]
rownames(ad_meta)=ad_meta$sample
ad_meta$sample=NULL
ad_taxa=ad_taxa[order(rownames(ad_taxa)),]
ad_meta=ad_meta[order(rownames(ad_meta)),]
rownames(ad_taxa)==rownames(ad_meta)
ad_meta$patient=as.factor(ad_meta$patient)

ad = foreach(i=1:ncol(ad_meta),.combine = rbind) %do%  {
  tmp.cov=ad_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=ad_taxa[rownames(ad_taxa) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.ad=adonis(tmp.taxa ~ tmp.cov[,1] ,  permutations = 1000, method = "euclidian") # or bray, or different methods
  
  cat(green(colnames(ad_meta)[i],"\n"))
  
  return.string=data.frame(Cov=colnames(ad_meta)[i],DF=tmp.ad$aov.tab[1,1],R2=tmp.ad$aov.tab[1,5],Pvalue=tmp.ad$aov.tab[1,6])
  
}
ad=ad[order(ad$R2,decreasing = F),]
rownames(ad)=ad$Cov
ad$Cov=NULL
order=factor(rownames(ad),levels = rownames(ad))
ad$Significant="No"
ad$Significant[ad$Pvalue<0.05]="Yes"
ggplot(data=ad, aes(x=order, y=R2,fill=Significant)) +
  geom_bar(stat="identity")+scale_fill_jama()+
  coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+xlab("")+guides(fill=F)
ggsave("OutputFigure/Adonis.pdf",width = 5,height = 2.5)

# composition plot
genus_plot=CompositionTable(data_bacteria,20)

tmp.id=as.character(genus_plot$ID[genus_plot$ID %in% c(data_meta$sample[data_meta$Group=="H"])])
genus_plot_H=AverageTable(genus_plot[genus_plot$ID %in% tmp.id,])
colnames(genus_plot_H)[1]="AverageAbundance"
genus_plot_H$Group="H"

tmp.id=as.character(genus_plot$ID[genus_plot$ID %in% c(data_meta$sample[data_meta$Group=="T"])])
genus_plot_T=AverageTable(genus_plot[genus_plot$ID %in% tmp.id,])
colnames(genus_plot_T)[1]="AverageAbundance"
genus_plot_T$Group="T"

tmp.id=as.character(genus_plot$ID[genus_plot$ID %in% c(data_meta$sample[data_meta$Group=="N"])])
genus_plot_N=AverageTable(genus_plot[genus_plot$ID %in% tmp.id,])
colnames(genus_plot_N)[1]="AverageAbundance"
genus_plot_N$Group="N"

genus_plot=rbind(genus_plot_H,genus_plot_T,genus_plot_N)

genus_plot_H=genus_plot_H[order(genus_plot_H$AverageAbundance,decreasing = T),]
Tax_s_order<-(genus_plot_H$Taxa)
palette=TaxaColor(Tax_s_order,level = "genus")
genus_plot$Taxa<-factor(genus_plot$Taxa,levels = rev(Tax_s_order))
ggplot(genus_plot, aes(x=Group, y=AverageAbundance,fill = Taxa))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
  ylab("Relative abundance")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+guides(fill=F)
ggsave("OutputFigure///Composition.genus.pdf",width = 4,height = 2)













































