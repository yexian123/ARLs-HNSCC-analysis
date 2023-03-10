#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)            
expFile="mRNA.txt"        
geneFile="gene.txt"       
setwd("C:\\4.ARGexp")    


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table(geneFile, header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="ARGexp.txt",sep="\t",quote=F,col.names=F)
