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

data_promoter=load("Input/启动子区域甲基化位点矩阵.RData")
data_enhancer=load("Input/增强子位点甲基化矩阵和位点甲基化注释基因.RData")
data_raw=read.table("Input/Candidates.promoter.matrix.txt",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
data_raw$ID=NULL
data_raw=as.data.frame(t(data_raw))
data_annot=read.table("Input/Site_420_name.addGeneName.from800w.txt",header = T,stringsAsFactors = F)
data_enhancer=as.data.frame(t(mWAS_meth))
data_promoter=data_raw[,colnames(data_raw) %in% data_annot$SampleID]

data_bacteria=read.table("Input/bacteria_matrix.txt",row.names = 1,sep = "\t",stringsAsFactors = F,check.names = F,header = T)
data_bacteria=as.data.frame(t(data_bacteria))
str(data_bacteria)
data_bacteria=data_bacteria[,colSums(data_bacteria!=0) > (nrow(data_bacteria) * 0.2) ]
data_bacteria=data_bacteria[,!colnames(data_bacteria) %in% c("Other","uncultured")]
data_bacteria=apply(data_bacteria,1,function(x){
  x=x/sum(x)
  return(x)
})
data_bacteria=as.data.frame(t(data_bacteria))

# re-arrange ID
data_meta=read.table("Input/Input.meta.txt",stringsAsFactors = F,header = T)
data_coupling=read.table("Input/clinical.txt",header = T,stringsAsFactors = F)
data_bacteria=data_bacteria[rownames(data_bacteria) %in% data_coupling$Tumor_Sample_Barcode,]
data_meta$sample=gsub("S","",data_meta$sample)

# step 1, associate bacteria with promoter

result_promoter = foreach(i=1:ncol(data_promoter),.combine = rbind) %do%  {
  
  cat(yellow(i,"\n"))
  tmp.probe=colnames(data_promoter)[i]
  tmp.data=merge(data_bacteria[,"Bacteroides",drop=F],data_coupling,by.x="row.names",by.y="Tumor_Sample_Barcode",all=F)
  colnames(tmp.data)[1]="Tumor_Sample_Barcode"
  rownames(tmp.data)=tmp.data$Tumor_Sample_Barcode
  tmp.data=merge(tmp.data,data_promoter[,tmp.probe,drop=F],by="row.names",all=F)
  tmp.data$Row.names=gsub("HCC","",tmp.data$Row.names)
  tmp.data=merge(tmp.data,data_meta,by.x="Row.names",by.y="sample",all=F)
  
  colnames(tmp.data)[3]="bacteria"
  colnames(tmp.data)[7]="probe"
  
  tmp.data$bacteria[tmp.data$bacteria==0]=NA
  tmp.data$bacteria=log2(tmp.data$bacteria)
  #tmp.data$bacteria=rmOutlier(tmp.data$bacteria)
  tmp.data=tmp.data[,c("Row.names","bacteria","probe","age","sex","BMI","smoking","drinking","Patient_ID")]
  tmp.data=tmp.data[tmp.data$Row.names %like% "T",]
  tmp.data=na.omit(tmp.data)
  Nonzero=as.numeric(length(tmp.data$probe[tmp.data$probe>0]))
  
  #tmp.data$probe[tmp.data$probe>0]=1
  tmp.data$probe=tmp.data$probe/100
  tmp.data$probe[tmp.data$probe==0]=0.001 # impute minimum
  tmp.data$probe[tmp.data$probe==1]=0.999 # impute maximum
  tmp.data$probe=log2(tmp.data$probe/(1-tmp.data$probe))
  
  if(length(unique(tmp.data$probe))>1){
    mod=lme(bacteria ~ probe + age+sex+BMI+smoking+drinking, random = ~ 1|Patient_ID,data = tmp.data)
    ano=anova(mod)
    mod=summary(mod)
    
    return.string=data.frame(Bacteria="Bacteroides",Probe=tmp.probe,beta=mod$tTable[,1][2],se=mod$tTable[,2][2],pvalue=ano$`p-value`[2],ProbeNonZero=Nonzero)
  }
  
}
result_promoter=merge(result_promoter,data_annot,by.x="Probe",by.y="SampleID",all=F)

# step 2, associate bacteria with enhancer

result_enhancer = foreach(i=1:ncol(data_enhancer),.combine = rbind) %do%  {
  
  cat(yellow(i,"\n"))
  tmp.probe=colnames(data_enhancer)[i]
  tmp.data=merge(data_bacteria[,"Bacteroides",drop=F],data_coupling,by.x="row.names",by.y="Tumor_Sample_Barcode",all=F)
  colnames(tmp.data)[1]="Tumor_Sample_Barcode"
  rownames(tmp.data)=tmp.data$Tumor_Sample_Barcode
  tmp.data=merge(tmp.data,data_enhancer[,tmp.probe,drop=F],by="row.names",all=F)
  tmp.data$Row.names=gsub("HCC","",tmp.data$Row.names)
  tmp.data=merge(tmp.data,data_meta,by.x="Row.names",by.y="sample",all=F)
  
  colnames(tmp.data)[3]="bacteria"
  colnames(tmp.data)[7]="probe"
  
  tmp.data$bacteria[tmp.data$bacteria==0]=NA
  tmp.data$bacteria=log2(tmp.data$bacteria)
  #tmp.data$bacteria=rmOutlier(tmp.data$bacteria)
  tmp.data=tmp.data[,c("Row.names","bacteria","probe","age","sex","BMI","smoking","drinking","Patient_ID")]
  tmp.data=tmp.data[tmp.data$Row.names %like% "T",]
  tmp.data=na.omit(tmp.data)
  Nonzero=as.numeric(length(tmp.data$probe[tmp.data$probe>0]))

  if(length(unique(tmp.data$probe))>1){
    mod=lme(bacteria ~ probe + age+sex+BMI+smoking+drinking, random = ~ 1|Patient_ID,data = tmp.data)
    ano=anova(mod)
    mod=summary(mod)
    
    return.string=data.frame(Bacteria="Bacteroides",Probe=tmp.probe,beta=mod$tTable[,1][2],se=mod$tTable[,2][2],pvalue=ano$`p-value`[2],ProbeNonZero=Nonzero)
  }
  
}
result_enhancer=merge(result_enhancer,nearGenes,by.x="Probe",by.y="ID",all=F)


























































