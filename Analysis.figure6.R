
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
rmOutlier=function(x){
  tmp.sd=sd(x,na.rm=T)
  tmp.rm=c(x[x<median(x,na.rm = T)-tmp.sd], x[x>median(x,na.rm = T)+tmp.sd])
  x[x %in% tmp.rm]=NA
  return(x)
}

# import data
data_bacteria=read.table("Input//bacteria_matrix.txt",row.names = 1,header = T,stringsAsFactors = F)
data_bacteria=as.data.frame(t(data_bacteria))
data_bacteria_log=data_bacteria
data_bacteria_log[data_bacteria_log==0]=NA
data_bacteria_log=log2(data_bacteria_log)
data_bacteria_log=data_bacteria_log[rownames(data_bacteria_log) %like% "T",]

data_methylation_candidate=read.table("Input/Gene420.from800w.df.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
data_methylation_candidate=as.data.frame(t(data_methylation_candidate))
data_methylation_candidate=data_methylation_candidate[,colSums(!is.na(data_methylation_candidate))>(0.9*nrow(data_methylation_candidate))]
data_methylation_candidate=data_methylation_candidate/100

data_methylation_candidate=data_methylation_candidate[rownames(data_methylation_candidate) %like% "T",]
data_methylation_candidate[data_methylation_candidate==0]=min(data_methylation_candidate[data_methylation_candidate!=0],na.rm =T)
data_methylation_candidate[data_methylation_candidate==1]=max(data_methylation_candidate[data_methylation_candidate!=1],na.rm =T)

data_methylation_candidate=log2(data_methylation_candidate/(1-data_methylation_candidate))
data_methylation_candidate=data_methylation_candidate[,-nearZeroVar(data_methylation_candidate)]

coupling_bacteria=read.table("Input/bacteria_sample_clinical.txt",header = T,stringsAsFactors = F)
coupling_methylation=read.table("Input/meth_sample_clinical.txt",header = T,stringsAsFactors = F)
coupling_bacteria$patient=paste("HCC",coupling_bacteria$patient,sep = "")

data_gene=read.table("Input/303_rpkm.txt")
data_gene=as.data.frame(t(data_gene))
data_gene=data_gene[,colSums(data_gene>0)>(nrow(data_gene)*0.8)]
data_gene[data_gene==0]=NA
data_gene_log=log2(data_gene) # log
protein_coding=read.table("Input/Protein.coding.list.ensemble.txt",header = T,stringsAsFactors = F)
data_gene_log=data_gene_log[,colnames(data_gene_log) %in% protein_coding$hgnc_symbol]

# CCA
tmp.gene=data_gene_log[rownames(data_gene_log) %in% rownames(data_bacteria_log),]
tmp.gene=tmp.gene[order(rownames(tmp.gene)),]
tmp.bacteria=data_bacteria_log[order(rownames(data_bacteria_log)),]
rownames(tmp.bacteria)==rownames(tmp.gene)
tmp.bacteria[is.na(tmp.bacteria)]=0
tmp.gene[is.na(tmp.gene)]=0
ccamodel_gene <- rda(tmp.gene~., tmp.bacteria)
plot(ccamodel_gene)

tmp.meth=data_methylation_log[rownames(data_methylation_log) %in% rownames(data_bacteria_log),]
tmp.meth=tmp.meth[order(rownames(tmp.meth)),]
tmp.bacteria=data_bacteria_log[rownames(data_bacteria_log) %in% rownames(tmp.meth),]
tmp.bacteria=tmp.bacteria[order(rownames(tmp.bacteria)),]
rownames(tmp.bacteria)==rownames(tmp.meth)
tmp.bacteria[is.na(tmp.bacteria)]=0
tmp.meth=tmp.meth[,-nearZeroVar(tmp.meth)]
ccamodel_meth <- rda(tmp.meth ~., tmp.bacteria)
plot(ccamodel_meth)

# adonis #### 
beta_diversity=vegdist(data_bacteria[rownames(data_bacteria) %like% "T",],method = "bray")
pcoa_analysis=as.data.frame(cmdscale(beta_diversity,k=10))

tmp.gene=data_gene_log[rownames(data_gene_log) %in% rownames(pcoa_analysis),]
tmp.gene=tmp.gene[order(rownames(tmp.gene)),]
tmp.bacteria=pcoa_analysis[order(rownames(pcoa_analysis)),]
rownames(tmp.bacteria)==rownames(tmp.gene)

adonis_gene=foreach(i=1:ncol(tmp.bacteria),.combine = rbind) %do%  {
  tmp.dim=colnames(tmp.bacteria)[i]
  
  tmp.result=adonis(formula = tmp.gene ~ tmp.bacteria[,tmp.dim],permutations = 100, method = "euclidean",na.rm = T)
  tmp.result=tmp.result$aov.tab
  
  cat(tmp.dim,"\n")
  return.string=data.frame(Dim=tmp.dim,R2=tmp.result$R2[1],P=tmp.result$`Pr(>F)`[1])
}
sum(adonis_gene$R2[adonis_gene$P<0.05])

tmp.meth=data_methylation_log[rownames(data_methylation_log) %in% rownames(pcoa_analysis),]
tmp.meth=tmp.meth[order(rownames(tmp.meth)),]
tmp.bacteria=pcoa_analysis[rownames(pcoa_analysis) %in% rownames(tmp.meth),]
tmp.bacteria=tmp.bacteria[order(rownames(tmp.bacteria)),]
rownames(tmp.bacteria)==rownames(tmp.meth)

adonis_meth=foreach(i=1:ncol(tmp.bacteria),.combine = rbind) %do%  {
  tmp.dim=colnames(tmp.bacteria)[i]
  
  tmp.result=adonis(formula = tmp.meth ~ tmp.bacteria[,tmp.dim],permutations = 100, method = "euclidean")
  tmp.result=tmp.result$aov.tab
  
  cat(tmp.dim,"\n")
  return.string=data.frame(Dim=tmp.dim,R2=tmp.result$R2[1],P=tmp.result$`Pr(>F)`[1])
}
sum(adonis_meth$R2[adonis_meth$P<0.05])

# gene-bacteria
# function of lmer
LmerTest = function(data_gene_log,tmp.gene,tmp.bacteria) {
  tmp.data=merge(data_bacteria_log[,tmp.bacteria,drop=F],coupling_bacteria,by.x="row.names",by.y="sample",all=F)
  tmp.data=merge(tmp.data,data_gene_log[,tmp.gene,drop=F],by.x="Row.names",by.y="row.names",all=F)
  colnames(tmp.data)=c("sample","bacteria","patient","gene")
  tmp.data=tmp.data[tmp.data$gene!="-Inf",]
  tmp.data=na.omit(tmp.data)
  
  mod = (lmer(gene~bacteria+(1|patient) +age + sex + BMI + smoking + drinking,data=tmp.data))
  mod.coef=as.data.frame(summary(mod)$coefficient)
  mod.sig=anova(mod)
  cat(green(tmp.bacteria,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,Bacteria=tmp.bacteria,Estimate=mod.coef$Estimate[2],SE=mod.coef$`Std. Error`[2],Pvalue=mod.sig$`Pr(>F)`)
}
tmp.bacteria="Bacteroides"
system.time(ind.res <- mclapply(mc.preschedule=F,mc.cores=8,
                                setNames(seq_len(ncol(data_gene_log)), dimnames(data_gene_log)[[2]]), LmerTest, data_gene_log=data_gene_log, 
                                tmp.bacteria=tmp.bacteria))
result_bacterioides <- rbindlist(ind.res,idcol="GeneSymbol")
result_bacterioides$FDR=p.adjust(result_bacterioides$Pvalue,method = "BH")
result_bacterioides_sig=result_bacterioides[result_bacterioides$Pvalue<0.05,]

result_bacterioides_enrich=read.table("OutputTable /result.Reactome.gene.down.txt",header = T,sep = "\t",stringsAsFactors = F,fill=T,quote = "",check.names = F)
result_bacterioides_enrich=result_bacterioides_enrich[order(result_bacterioides_enrich$`Entities FDR`,decreasing = F),]
result_bacterioides_enrich=result_bacterioides_enrich[1:10,]
result_bacterioides_enrich$`Pathway name`=factor(result_bacterioides_enrich$`Pathway name`,levels = rev(result_bacterioides_enrich$`Pathway name`))
result_bacterioides_enrich$`Entities FDR`=-log10(result_bacterioides_enrich$`Entities FDR`)
ggplot(result_bacterioides_enrich, aes(x = `Entities FDR`, y = `Pathway name`, size = `Entities found`)) +
  geom_point(color = "#AA4499", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave("OutputFigure/Gene.pathway.pdf",width = 7,height = 3)

# deconvolution
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

result_epic=foreach(i=1:ncol(data_bacteria),.combine = rbind) %do%  {
  tmp.bacteria=colnames(data_bacteria)[i]
  
  tmp.result=foreach(j=1:ncol(data_epic),.combine = rbind) %do%  {
    tmp.cell=colnames(data_epic)[j]
    tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],coupling_bacteria,by.x="row.names",by.y="sample",all=F)
    tmp.data=merge(tmp.data,data_epic[,tmp.cell,drop=F],by.x="Row.names",by.y="row.names",all=F)
    tmp.data=na.omit(tmp.data)
    colnames(tmp.data)=c("sample","bacteria","patient","cell")
    tmp.data$cell=log2(tmp.data$cell)
    tmp.data$cell=rmOutlier(tmp.data$cell)
    tmp.data=tmp.data[tmp.data$cell!="-Inf",]
    
    mod = (lmer(cell~bacteria+(1|patient) +age + sex + BMI + smoking + drinking,data=tmp.data))
    mod.coef=as.data.frame(summary(mod)$coefficient)
    mod.sig=anova(mod)
    cat(green(tmp.bacteria,"+++",tmp.cell,"\n"))
    
    return.string=data.frame(Cell=tmp.cell,Bacteria=tmp.bacteria,Estimate=mod.coef$Estimate[2],SE=mod.coef$`Std. Error`[2],Pvalue=mod.sig$`Pr(>F)`)
  }
  return.string=tmp.result
}
write.table(result_epic,"OutputTable /Bacteroides.epic.cell.txt",sep = "\t",row.names = F,quote = F)

result_ciber=foreach(i=1:ncol(data_bacteria),.combine = rbind) %do%  {
  tmp.bacteria=colnames(data_bacteria)[i]
  
  tmp.result=foreach(j=1:ncol(data_ciber),.combine = rbind) %do%  {
    tmp.cell=colnames(data_ciber)[j]
    tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],coupling_bacteria,by.x="row.names",by.y="sample",all=F)
    tmp.data=merge(tmp.data,data_ciber[,tmp.cell,drop=F],by.x="Row.names",by.y="row.names",all=F)
    tmp.data=na.omit(tmp.data)
    colnames(tmp.data)=c("sample","bacteria","patient","cell")
    tmp.data$cell=log2(tmp.data$cell)
    tmp.data$cell=rmOutlier(tmp.data$cell)
    tmp.data=tmp.data[tmp.data$cell!="-Inf",]
    tmp.data=na.omit(tmp.data)
    
    if(sum(tmp.data$cell)!=0){
      mod = (lmer(cell~bacteria+(1|patient) +age + sex + BMI + smoking + drinking,data=tmp.data))
      mod.coef=as.data.frame(summary(mod)$coefficient)
      mod.sig=anova(mod)
      cat(green(tmp.bacteria,"+++",tmp.cell,"\n"))
      
      return.string=data.frame(Cell=tmp.cell,Bacteria=tmp.bacteria,Estimate=mod.coef$Estimate[2],SE=mod.coef$`Std. Error`[2],Pvalue=mod.sig$`Pr(>F)`)
    }
  }
  return.string=tmp.result
}
write.table(result_ciber,"OutputTable /Bacteroides.ciber.cell.txt",sep = "\t",row.names = F,quote = F)

result_xcell=foreach(i=1:ncol(data_bacteria),.combine = rbind) %do%  {
  tmp.bacteria=colnames(data_bacteria)[i]
  
  tmp.result=foreach(j=1:ncol(data_xcell),.combine = rbind) %do%  {
    tmp.cell=colnames(data_xcell)[j]
    tmp.data=merge(data_bacteria[,tmp.bacteria,drop=F],coupling_bacteria,by.x="row.names",by.y="sample",all=F)
    tmp.data=merge(tmp.data,data_xcell[,tmp.cell,drop=F],by.x="Row.names",by.y="row.names",all=F)
    tmp.data=na.omit(tmp.data)
    colnames(tmp.data)=c("sample","bacteria","patient","cell")
    tmp.data=tmp.data[tmp.data$cell!=0,]
    tmp.data$cell=log2(tmp.data$cell)
    tmp.data$cell=rmOutlier(tmp.data$cell)
    tmp.data=tmp.data[tmp.data$cell!="-Inf",]
    tmp.data=na.omit(tmp.data)
    
    mod = (lmer(cell~bacteria+(1|patient) +age + sex + BMI + smoking + drinking,data=tmp.data))
    mod.coef=as.data.frame(summary(mod)$coefficient)
    mod.sig=anova(mod)
    cat(green(tmp.bacteria,"+++",tmp.cell,"\n"))
    
    return.string=data.frame(Cell=tmp.cell,Bacteria=tmp.bacteria,Estimate=mod.coef$Estimate[2],SE=mod.coef$`Std. Error`[2],Pvalue=mod.sig$`Pr(>F)`)
  }
  return.string=tmp.result
}
write.table(result_xcell,"OutputTable /Bacteroides.ciber.cell.txt",sep = "\t",row.names = F,quote = F)


library(dplyr)
library(stringr)
load('Input/Rarefy_phyloseq.RData')
load('Input/细胞解卷积和转录组.RData')
library(IOBR)
clincal_final <- clincal_final[c(7,8)]
clincal_final$sample=str_remove(clincal_final$sample,"S")
clincal_final <- clincal_final %>% column_to_rownames('sample')
rownames(clincal_final)=paste("HCC",rownames(clincal_final),sep="")
sig <- signature_collection_citation[!duplicated(signature_collection_citation$Journal),]
eset_stad<-count2tpm(countMat = reads, idType = "SYMBOL", org="hsa", source = "local" )


sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_collection,
                             method          = "zscore",
                             mini_gene_count = 2)

tmp.sig <- sig_tme %>% column_to_rownames('ID')
data_taxa=data_bacteria_log
tmp.sig=tmp.sig[rownames(data_taxa),]

library(nlme)
result_mixedLm = foreach(i=1:ncol(tmp.sig),.combine = rbind) %do%  {
  tmp.cell=colnames(tmp.sig)[i]
  tmp.bacteria=data_taxa[,"Bacteroides",drop=F]
  tmp.data=merge(tmp.bacteria,tmp.sig[,tmp.cell,drop=F],by="row.names",all=F)
  tmp.data=merge(tmp.data,clincal_final,by.x="Row.names",by.y="row.names",all=F)
  colnames(tmp.data)=c("ID","Bacteria","Cell","Subject")
  
  tmp.data=tmp.data[tmp.data$Bacteria!=-Inf & tmp.data$Cell!=-Inf,]
  tmp.data=na.omit(tmp.data)
  
  mod=(lme(Cell ~ Bacteria, tmp.data, random = ~ 1|Subject))
  #mod=lm(Cell ~ Bacteria, tmp.data,)
  ano=anova(mod)
  mod=summary(mod)
  cat(yellow(i, "+++++", tmp.cell,"\n"))
  
  return.string=data.frame(cell=tmp.cell,Beta=mod$coefficients$fixed[2],Pvalue=ano$`p-value`[2])
  #return.string=data.frame(cell=tmp.cell,Beta=mod$coefficients[2],Pvalue=mod$coefficients[8])
}
result_mixedLm$FDR=p.adjust(result_mixedLm$Pvalue,method = "BH")
sig_up=as.character(result_mixedLm$cell[result_mixedLm$Pvalue<0.05 & result_mixedLm$Beta>0])
sig_down=as.character(result_mixedLm$cell[result_mixedLm$Pvalue<0.05 & result_mixedLm$Beta<0])

draw_df=cbind(tmp.sig,tmp.bacteria)

#IFG
ifg_data=draw_df %>% dplyr::select(Bacteroides,IFNG_signature_Ayers_et_al)
ifg_data=ifg_data[ifg_data$Bacteroides!=-Inf & ifg_data$IFNG_signature_Ayers_et_al!=-Inf,]
ggplot(data = ifg_data,aes(x=Bacteroides,y=IFNG_signature_Ayers_et_al))+geom_point(size=3,alpha=0.8)+
  geom_smooth(method = 'lm')+stat_cor(method = "spearman",label.x = -3)+theme_bw()

ggplot(ifg_data, aes(Bacteroides, IFNG_signature_Ayers_et_al)) +
  geom_point(shape = 21, color = "black", size = 2,fill="#AA4499")+
  geom_smooth(method = lm,color="grey")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+guides(color=F)

ggsave("OutputFigure/INFY.pdf",width = 3,height = 3)



















































