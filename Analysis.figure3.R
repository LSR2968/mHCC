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
library(tidyverse)
library(ggnewscale)
library(Hmisc)

load("Input/WES矩阵和临床信息表.RData")
## ------绘制这些患者的突变树和热图-------

library(MesKit)
library(phangorn)
## --------建立菌的进化树--------
file_all$patient=str_remove_all(file_all$Tumor_Sample_Barcode,"[TM].*")
file_all=file_all[file_all$Tumor_Sample_Barcode != "HCC817M",]
file_all2 <- file_all %>% dplyr::select(Tumor_Seq_Allele2,Chromosome,Start_Position,Reference_Allele,Tumor_Sample_Barcode,patient)
file_all2 <- unite(file_all2,'gene',c('Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2'),sep = "&")

wes_mat <-split(file_all2,file_all2$patient)

byNJ <- function(mut_dat){
  matTree <- nj(dist.gene(mut_dat))
  return(matTree)
}

wes_mat <-lapply(wes_mat,function(x){x=x %>% 
  reshape2::dcast(Tumor_Sample_Barcode~gene,valvue.vars=,funaggregate=toString) %>% column_to_rownames('Tumor_Sample_Barcode');
x =ifelse(is.na(x),0,1) %>% as.data.frame();})

wes_mat <- lapply(wes_mat,as.data.frame)
for (i in names(wes_mat)){
  wes_mat[[i]]=wes_mat[[i]] %>%   t() %>% as.data.frame() %>%  dplyr::mutate(!!paste0(i,'L'):=0) %>% t() %>% as.data.frame()
}
wes_tree_list <-   lapply(wes_mat,
                          function(x){y=byNJ(as.matrix(x))})

#----------导入菌群丰度矩阵---
bac_mat <- read.csv('Input/筛选过后的菌属丰度表_26个的.csv',row.names = 1,check.names = F)
colnames(bac_mat)=str_remove(colnames(bac_mat),"S")
colnames(bac_mat)=paste0('HCC',colnames(bac_mat))
bac_mat=t(bac_mat) %>% as.data.frame()
bac_mat=bac_mat[rownames(bac_mat) != "HCC817M",]

bac_mat_clr <- zCompositions::cmultRepl(bac_mat, method="CZM", label=0)
bac_mat_clr = compositions::clr(bac_mat_clr)
bac_mat_clr=as.data.frame(bac_mat_clr)

# whether to use clr transformed or 0/1 matrix
#bac_mat=bac_mat_clr

# 建立fitness函数----
evaluate_dist=function(indices=c()){
  a=try({
    group=str_remove(rownames(bac_mat),"[LTM].*")
    bac_list=split(bac_mat[,indices==1],group)
    bac_tree_list <- lapply(bac_list,
                            function(x){y=nj(vegdist(x,method = 'euclidian'))})
    mean_v=c()
    for (i in names(bac_list)){
      mean_v=c(mean_v,RF.dist(bac_tree_list[[i]],wes_tree_list[[i]]))
    }
  },silent = T)
  if ('try-error' %in% class(a)){
    return(-100000)
  }else{return(-mean(mean_v))}
}


# 运行GA算法----
library(GA)
woppa <- ga(
  type='binary',
  fitness = evaluate_dist,
  seed = 55,
  popSize = 10,
  nBits = ncol(bac_mat),
  parallel=10)

res <-   woppa@solution  %>% as.data.frame()
colnames(res)=colnames(bac_mat)
res <- res %>% dplyr::mutate(sum=rowSums(.)) %>% dplyr::arrange(sum)
colnames(res)[res[1,]==1]

group=str_remove(rownames(bac_mat),"[LT].*")

bac_list=split(bac_mat[,colnames(res)[res[1,]==1]],group)

bac_tree_list <- lapply(bac_list,
                        function(x){y=njs(vegdist(x,method = 'jaccard'))})

for (i in names(bac_tree_list)){
  bac_tree_list[[i]]=root(bac_tree_list[[i]],paste0(i,'L'),resolve.root = T)
}
mean_v=c()
for (i in names(bac_list)){
  mean_v=c(mean_v,RF.dist(bac_tree_list[[i]],wes_tree_list[[i]]))
}
mean_v
mean(mean_v)
table(mean_v)[1]/length(mean_v)
names(bac_list)[mean_v==0]
loop_real=data.frame(NumberPermutate=0,Bacteria=paste(c(colnames(res)[res[1,]==1]), collapse = " "),DistanceMean=mean(mean_v),ProP=table(mean_v)[1]/length(mean_v))

# permutation for significance
loop_permute = foreach(n=1:100,.combine = rbind) %do%  {
  set.seed(n)
  cat(yellow(n),"\n")
  
  tmp.bac_mat=bac_mat[sample(nrow(bac_mat)),]
  tmp.bac_mat=tmp.bac_mat[,sample(ncol(tmp.bac_mat))]
  colnames(tmp.bac_mat)=colnames(bac_mat)
  rownames(tmp.bac_mat)=rownames(bac_mat)
  
  evaluate_dist_permute=function(indices=c()){
    a=try({
      group=str_remove(rownames(tmp.bac_mat),"[LTM].*")
      bac_list=split(tmp.bac_mat[,indices==1],group)
      bac_tree_list <- lapply(bac_list,
                              function(x){y=nj(vegdist(x,method = 'jaccard'))})
      mean_v=c()
      for (i in names(bac_list)){
        mean_v=c(mean_v,RF.dist(bac_tree_list[[i]],wes_tree_list[[1]])) # random select WES tree, repeat 100 times, here the first tree for example
      }
    },silent = T)
    if ('try-error' %in% class(a)){
      return(-100000)
    }else{return(-mean(mean_v))}
  }
  woppa <- ga(
    type='binary',
    fitness = evaluate_dist_permute,
    seed = 1111,
    popSize = 100,
    nBits = ncol(tmp.bac_mat),
    parallel=10)
  
  res <-   woppa@solution  %>% as.data.frame()
  colnames(res)=colnames(tmp.bac_mat)
  res <- res %>% dplyr::mutate(sum=rowSums(.)) %>% dplyr::arrange(sum)
  colnames(res)[res[1,]==1]
  
  group=str_remove(rownames(tmp.bac_mat),"[LT].*")
  
  bac_list=split(tmp.bac_mat[,colnames(res)[res[1,]==1]],group)
  
  bac_tree_list <- lapply(bac_list,
                          function(x){y=njs(vegdist(x,method = 'jaccard'))})
  
  for (i in names(bac_tree_list)){
    bac_tree_list[[i]]=root(bac_tree_list[[i]],paste0(i,'L'),resolve.root = T)
  }
  mean_v=c()
  for (i in names(bac_list)){
    mean_v=c(mean_v,RF.dist(bac_tree_list[[i]],wes_tree_list[[i]]))
  }
  mean_v
  names(bac_list)[mean_v==0]
  
  return.string=data.frame(NumberPermutate=n,Bacteria=paste(c(colnames(res)[res[1,]==1]), collapse = " "),DistanceMean=mean(mean_v),ProP=table(mean_v)[1]/length(mean_v))
  #return.string=data.frame(NumberPermutate=n,Patient_dis=mean_v)
  }

#loop_permute=read.table("OutputTable /Loop_permutation.binary.txt",header = T,sep="\t")
loop_plot=rbind(loop_real,loop_permute)
loop_plot$NumberPermutate=as.factor(loop_plot$NumberPermutate)
length(loop_plot$ProP[loop_plot$ProP<=loop_real$ProP])
#write.table(loop_plot,file = "OutputTable //Loop_permutation.binary.txt",sep = "\t",quote = F,row.names = F)

ggplot(loop_permute, aes(x=ProP)) +
  geom_density(color="#999999",bw=0.02)+geom_vline(xintercept = loop_real$ProP,linetype="dotted")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+guides(fill=F)
ggsave("OutputFigure/Supple.X.pdf",width = 3,height = 3)

# repeated test for Bacteroides
loop_repeat=data.frame(Bacteria=colnames(bac_mat),RepeatedNumber=NA,CountNumber=0)
for(n in 1:100){
  set.seed(n)
  cat(yellow(n),"\n")
  
  woppa <- ga(
    type='binary',
    fitness = evaluate_dist,
    seed = n,
    popSize = 5,
    nBits = ncol(bac_mat),
    parallel=10)
  
  res <-   woppa@solution  %>% as.data.frame()
  colnames(res)=colnames(bac_mat)
  res <- res %>% dplyr::mutate(sum=rowSums(.)) %>% dplyr::arrange(sum)
  tmp.selected=colnames(res)[res[1,]==1]
  
  for(i in 1:nrow(loop_repeat)){
    tmp.bacteria=loop_repeat$Bacteria[i]
    if(tmp.bacteria %in% tmp.selected){
      loop_repeat$CountNumber[loop_repeat$Bacteria==tmp.bacteria]=loop_repeat$CountNumber[loop_repeat$Bacteria==tmp.bacteria]+1
    }
  }
}
write.table(loop_repeat,file = "OutputTable /Loop_repeat.binary.txt",row.names = F,quote = F,sep = "\t")

loop_repeat=data.frame(Bacteria=colnames(bac_mat),RepeatedNumber=NA,CountNumber=0)
for(n in 1:100){
  set.seed(n)
  cat(yellow(n),"\n")
  
  woppa <- ga(
    type='binary',
    fitness = evaluate_dist,
    seed = n,
    popSize = 5,
    nBits = ncol(bac_mat),
    parallel=10)
  
  res <-   woppa@solution  %>% as.data.frame()
  colnames(res)=colnames(bac_mat)
  res <- res %>% dplyr::mutate(sum=rowSums(.)) %>% dplyr::arrange(sum)
  tmp.selected=colnames(res)[res[1,]==1]
  
  for(i in 1:nrow(loop_repeat)){
    tmp.bacteria=loop_repeat$Bacteria[i]
    if(tmp.bacteria %in% tmp.selected){
      loop_repeat$CountNumber[loop_repeat$Bacteria==tmp.bacteria]=loop_repeat$CountNumber[loop_repeat$Bacteria==tmp.bacteria]+1
    }
  }
}
write.table(loop_repeat,file = "OutputTable /Loop_repeat.quantitiative.txt",row.names = F,quote = F,sep = "\t")

loop_repeat=read.table("OutputTable /Loop_repeat.binary.txt",header = T,stringsAsFactors = F)
loop_repeat=loop_repeat[loop_repeat$Bacteria %in% c("Acinetobacter","Bacteroides","Chryseobacterium","Corynebacterium",
                                                    "Enterococcus","Lawsonella","Muribaculaceae","Pseudomonas","Sphingomonas","Streptococcus"),]
#loop_repeat=loop_repeat[order(loop_repeat$CountNumber,decreasing = F),]
loop_repeat$Bacteria=factor(loop_repeat$Bacteria,levels = rev(c("Bacteroides","Lawsonella","Corynebacterium","Streptococcus",
                                                            "Chryseobacterium","Acinetobacter","Muribaculaceae","Sphingomonas","Pseudomonas","Enterococcus")))
ggplot(data=loop_repeat, aes(x=Bacteria, y=CountNumber)) +
  geom_bar(stat="identity",fill="#E7B800")+
  coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+xlab("")+guides(fill=F)
ggsave("OutputFigure/Supple.Y1.pdf",width = 4,height = 3)

# correlation bewteen 10 bacteria selected by GA
bac_mat_binary=bac_mat
bac_mat_binary[bac_mat_binary>0]=1

bac_mat_quantita=bac_mat
bac_mat_quantita[bac_mat_quantita==0]=NA
bac_mat_quantita=log2(bac_mat_quantita)

bac_mat_quantita.r=as.data.frame(cor(bac_mat_quantita,method = "spearman",use = "pairwise.complete.obs"))
bac_mat_binary.r=as.data.frame(cor(bac_mat_binary,method = "spearman"))
bac_mat_quantita.p=as.data.frame(rcorr(as.matrix(bac_mat_quantita),type = "spearman")$P)
bac_mat_binary.p=as.data.frame(rcorr(as.matrix(bac_mat_binary),type = "spearman")$P)

tmp.bacteria=c("Bacteroides","Lawsonella","Corynebacterium","Streptococcus","Chryseobacterium","Acinetobacter","Muribaculaceae","Sphingomonas","Pseudomonas","Enterococcus")
bac_mat_binary.r=bac_mat_binary.r[colnames(bac_mat_binary.r) %in% tmp.bacteria,rownames(bac_mat_binary.r) %in% tmp.bacteria]
bac_mat_quantita.r=bac_mat_quantita.r[colnames(bac_mat_quantita.r) %in% tmp.bacteria,rownames(bac_mat_quantita.r) %in% tmp.bacteria]
bac_mat_quantita.p=bac_mat_quantita.p[colnames(bac_mat_quantita.p) %in% tmp.bacteria,rownames(bac_mat_quantita.p) %in% tmp.bacteria]
bac_mat_binary.p=bac_mat_binary.p[colnames(bac_mat_binary.p) %in% tmp.bacteria,rownames(bac_mat_binary.p) %in% tmp.bacteria]

get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

get_lower_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}

upper_tri <- get_upper_tri(bac_mat_quantita.r)
lower_tri <- get_lower_tri(bac_mat_binary.r)
library(reshape2)
meltUp=gather(upper_tri, v1, cor)
meltUp$v2=rep(rownames(upper_tri),times=10) 

meltLow=gather(lower_tri, v1, cor)
meltLow$v2=rep(rownames(lower_tri),times=10) 
meltUp=na.omit(meltUp)
meltLow=na.omit(meltLow)
meltUp[meltUp<0]=0
meltLow[meltLow<0]=0

ggplot() +
  geom_tile(data=meltUp, aes(x=v1, y=v2, fill=cor))+
  scale_fill_gradient(low = "white", high = "#009988",
                      limit = c(0,1), name = "Correlation")+
  new_scale_fill() +
  geom_tile(data=meltLow, aes(x=v1, y=v2, fill=cor))+
    scale_fill_gradient(low = "white", high = "#E7B800",
                        limit = c(0,1), name = "Correlation")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  coord_fixed()+scale_x_discrete(guide = guide_axis(angle = 45)) +xlab("")+ylab("")
ggsave("OutputFigure/10.bacteria.cor.pdf",width = 5,height = 5)
























































