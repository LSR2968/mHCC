library(data.table)
library(vegan)
library("lmerTest")
library("languageR")
library(caret)
library(glmnet)
library(parallel)
library(crayon)
library(factoextra)
library(foreach)
library(nlme)
require(lme4)
library("ggplot2")
library("ggdendro")
library(randomcoloR)
library(data.table)

colorlist <- c("#30ACEC", '#1E9D89','#FB9A99','#FECD90','#B39DDB',
               "#00B9AF",'#E31A1C','#FF7F00','#2CA02C','#6A3D9A',
               '#00639C','#B9DE28','#FFEA46','#FABFD2','#DD53B4',
               '#1A66FF','#86FF86','#FFFF99','#B54923','#9900CC',
               '#C4D8F3','#D4EDCB','#F36D60','#F5B355','#A65CA4')
genus <- c("Other","Vibrio","Stenotrophomonas","Pseudomonas","Corynebacterium",
           "Methylobacterium-Methylorubrum","Acinetobacter" ,"Bacillus","Escherichia-Shigella","Prevotella",
           "Cutibacterium" ,"Streptococcus" ,"Enterococcus","Staphylococcus","Sphingomonas",
           "Muribaculaceae" ,"Bacteroides" ,"Lactobacillus" ,"Fusobacterium" ,"Bifidobacterium",
           "Chryseobacterium","Enhydrobacter" ,"Lawsonella"  ,"Paracoccus", "Brevundimonas")
names(colorlist) <- genus

data_bacteria=read.table("Input//bacteria_matrix.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
data_bacteria=as.data.frame(t(data_bacteria))
data_bacteria=data_bacteria[,colSums(data_bacteria>0)>(nrow(data_bacteria)*0.2)]
data_bacteria=data_bacteria[,!colnames(data_bacteria) %in% c("Other","uncultured")]
data_bacteria=apply(data_bacteria,1,function(x){
  x=x/sum(x)
  return(x)
})
data_bacteria=as.data.frame(t(data_bacteria))
data_bacteria_log=data_bacteria
data_bacteria_log[data_bacteria_log==0]=NA
data_bacteria_log=log2(data_bacteria_log)
data_bacteria_log=data_bacteria_log[rownames(data_bacteria_log) %like% "T",]

data_cell=load("Input/细胞解卷积和转录组.RData")

data_epic=as.data.frame(epic)
data_ciber=as.data.frame(cibersort)
data_xcell=as.data.frame(xcell)

rownames(data_epic)=data_epic$ID
data_epic$ID=NULL
rownames(data_ciber)=data_ciber$ID
data_ciber$ID=NULL
rownames(data_xcell)=data_xcell$ID
data_xcell$ID=NULL

data_epic=data_epic[,colSums(data_epic!=0)>(nrow(data_epic)*0.2)]
data_ciber=data_ciber[,colSums(data_ciber!=0)>(nrow(data_ciber)*0.2)]
data_xcell=data_xcell[,colSums(data_xcell!=0)>(nrow(data_xcell)*0.2)]
data_ciber=data_ciber[,1:17]

data_meta=read.table("Input/Input.meta.txt",stringsAsFactors = F,header = T)
data_meta$sample=gsub("S","",data_meta$sample)

data_gene=read.table("Input/303_rpkm.txt")
data_gene=as.data.frame(t(data_gene))
data_gene=data_gene[,colSums(data_gene>0)>(nrow(data_gene)*0.8)]
data_gene[data_gene==0]=NA
data_gene_log=log2(data_gene)
protein_coding=read.table("Input/Protein.coding.list.ensemble.txt",header = T,stringsAsFactors = F)
data_gene_log=data_gene_log[,colnames(data_gene_log) %in% protein_coding$hgnc_symbol]

data_gene_log=data_gene_log[rownames(data_gene_log) %like% "T",]
data_ciber=data_ciber[rownames(data_ciber) %in% rownames(data_gene_log),]
data_gene_log=data_gene_log[rownames(data_gene_log) %in% rownames(data_ciber),]
distance <- dist(scale(data_gene_log), method = "euclidean")
hc <- hclust(distance, method = "ward.D2")
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type = "rectangle")
mm=dend_data$segments
nn=dend_data$labels
nn$label=gsub("HCC","",nn$label)
nn=merge(nn,data_meta,by.x="label",by.y = "sample",all=F)
nn$dot=1
set.seed(1)
palette <- distinctColorPalette(58)
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_point(data=nn,aes(x,y,color=patient),
             size = 2, shape = 16) +scale_color_manual(values = palette)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+guides(color=F)
ggsave("OutputFigure/Tree.gene.plot.pdf",width = 15,height = 3)

data_cell=data_ciber
data_cell=stack(data_cell)
data_cell$patient=rep(rownames(data_ciber),n=22)
data_cell$patient=gsub("HCC","",data_cell$patient)
data_cell$patient=factor(data_cell$patient,levels = nn$label)
set.seed(1)
ggplot(data_cell, aes(x = patient, y = values, fill = ind)) +  # Create stacked bar chart
  geom_bar(stat = "identity")+scale_fill_manual(values = palette)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+guides(fill=F)
ggsave("OutputFigure/bar.gene.plot.pdf",width = 15,height = 3)

# assess variation of genes and cells within tumor, compared with GAPDH/ACTB
tmp.standard=c()
tmp.data=data_gene_log[,"ACTB",drop=F]
rownames(tmp.data)=gsub("HCC","",rownames(tmp.data))
tmp.gapdh=merge(tmp.data,data_meta,by.x="row.names",by.y="sample",all=F)
colnames(tmp.gapdh)[2]="GAPDH"
for(n in unique(tmp.gapdh$patient)){
  tmp.patient=tmp.gapdh[tmp.gapdh$patient==n,]
  tmp.mm=sd(tmp.patient$GAPDH) / mean(tmp.patient$GAPDH) 
  tmp.standard=append(tmp.standard,tmp.mm)
}
tmp.standard=mean(tmp.standard)

result_cv_gene = foreach(i=1:ncol(data_gene_log),.combine = cbind) %do%  {
  cat(i,"\n")
  
  tmp.gene=colnames(data_gene_log)[i]
  tmp.data=data_gene_log[,tmp.gene,drop=F]
  rownames(tmp.data)=gsub("HCC","",rownames(tmp.data))
  tmp.data=merge(tmp.data,data_meta,by.x="row.names",by.y="sample",all=F)
  colnames(tmp.data)[2]="Gene"
  
  tmp.cv=c()
  for(n in unique(tmp.data$patient)){
    tmp.patient=tmp.data[tmp.data$patient==n,1:3]
    tmp.patient=merge(tmp.patient,tmp.gapdh,by="Row.names",all=F)
    tmp.patient=tmp.patient[,c("GAPDH","Gene")]
    tmp.patient=stack(tmp.patient)
    tmp.patient=na.omit(tmp.patient)
    
    if(table(tmp.patient$ind)[2]>1){
      tmp.mod=var.test(values ~ ind, tmp.patient, 
                       alternative = "two.sided")
      tmp.cv=append(tmp.cv,tmp.mod$p.value)
    }else{tmp.cv=append(tmp.cv,NA)}

  }
  tmp.result=data.frame(tmp.cv,row.names = unique(tmp.data$patient))
  colnames(tmp.result)=tmp.gene
  
  return.string=tmp.result
  
}
mm=apply(result_cv_gene,2,function(x){
  tmp.count=length(na.omit(x[x<0.05]))
  return(tmp.count/length(x[!is.na(x)]))
})
mm=data.frame(Gene=colnames(data_gene_log),Aver_prop=mm)

result_cv_cell = foreach(i=1:ncol(data_xcell),.combine = cbind) %do%  {
  cat(i,"\n")
  
  tmp.gene=colnames(data_xcell)[i]
  tmp.data=data_xcell[,tmp.gene,drop=F]
  rownames(tmp.data)=gsub("HCC","",rownames(tmp.data))
  tmp.data=merge(tmp.data,data_meta,by.x="row.names",by.y="sample",all=F)
  colnames(tmp.data)[2]="Gene"
  tmp.data$Gene[tmp.data$Gene==0]=NA
  tmp.data$Gene=log2(tmp.data$Gene)
  
  tmp.cv=c()
  for(n in unique(tmp.data$patient)){
    tmp.patient=tmp.data[tmp.data$patient==n,1:3]
    tmp.patient=merge(tmp.patient,tmp.gapdh,by="Row.names",all=F)
    
    if(length(table(tmp.patient$Gene))>1){
      tmp.mod=sd(tmp.patient$Gene)/mean(tmp.patient$Gene)
      tmp.cv=append(tmp.cv,tmp.mod)
    }else{tmp.cv=append(tmp.cv,NA)}
    
  }
  tmp.result=data.frame(tmp.cv,row.names = unique(tmp.data$patient))
  colnames(tmp.result)=tmp.gene
  
  return.string=tmp.result
  
}
mm=apply(result_cv_cell,2,function(x){
  x=mean(x[!is.na(x)])
  return(x)
})
mm=data.frame(Gene=colnames(data_xcell),Aver_prop=mm)

calculate.dist=function(ref,alt,data){
  tmp.dist=c()
  tmp.ref=data[data$Row.names %like% ref,]
  tmp.alt=data[data$Row.names %like% alt,]
  for(i in tmp.alt$Row.names){
    tmp.mod=tmp.alt$Gene[tmp.alt$Row.names==i]-tmp.ref$Gene
    tmp.dist=append(tmp.dist,tmp.mod)
  }
  return(sd(tmp.dist)/mean(tmp.dist))
}
result_cv_bacteria = foreach(i=1:ncol(data_bacteria),.combine = cbind) %do%  {
  cat(i,"\n")
  
  tmp.gene=colnames(data_bacteria)[i]
  tmp.data=data_bacteria[,tmp.gene,drop=F]
  rownames(tmp.data)=gsub("HCC","",rownames(tmp.data))
  tmp.data=merge(tmp.data,data_meta,by.x="row.names",by.y="sample",all=F)
  colnames(tmp.data)[2]="Gene"
  
  tmp.cv=c()
  for(n in unique(tmp.data$patient)){
    tmp.patient=tmp.data[tmp.data$patient==n,1:3]
    
    if(length(table(tmp.patient$Gene))>2){
      tmp.patient=tmp.patient[tmp.patient$Row.names %like% "T",]
      tmp.mod=sd(tmp.patient$Gene)/mean(tmp.patient$Gene)
      tmp.cv=append(tmp.cv,tmp.mod)
    }else{tmp.cv=append(tmp.cv,NA)}
    
  }
  tmp.result=data.frame(tmp.cv,row.names = unique(tmp.data$patient))
  colnames(tmp.result)=tmp.gene
  
  return.string=tmp.result
  
}
mm=apply(result_cv_bacteria,2,function(x){
  x=median(x[!is.na(x)])
  return(x)
})
mm=data.frame(Gene=colnames(data_bacteria_log),Aver_prop=mm)

result_cv_bacteria=read.table("OutputTable /CV.bacteria.result.txt",sep = "\t",header = T,stringsAsFactors = F)
mm=read.table("OutputTable /CV.bacteria.summary.txt",header = T,stringsAsFactors = F)
mm=mm[order(mm$Aver_prop,decreasing = T),]
mm=mm[1:10,] # variation top 10
result_cv_bacteria=result_cv_bacteria[,colnames(result_cv_bacteria) %in% mm$Gene]
result_cv_bacteria=stack(result_cv_bacteria)

result_cv_bacteria$ind=factor(result_cv_bacteria$ind,levels = rev(mm$Gene))
ggplot(result_cv_bacteria,                          # Change colors of jitter points
       aes(x = ind,
           y = values,
           fill = ind)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 3,
                                             dodge.width = 0.7),
             aes(fill = ind),
             pch = 21)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+coord_flip()+scale_fill_manual(values = colorlist)+guides(fill=F)
































