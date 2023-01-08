#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.5            
pvalueFilter=0.001       
setwd("C:\\5.ARGlncExp")     


rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=read.table("ARGexp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
ARG=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
ARG=avereps(ARG)
ARG=ARG[rowMeans(ARG)>0.1,]


group=sapply(strsplit(colnames(ARG),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
ARG=ARG[,group==0]


outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(ARG)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(ARG[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(ARG=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(ARG=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}


write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)

lncNode=data.frame(Node=unique(as.vector(outTab[,"lncRNA"])), Type="lncRNA")
mrnaNode=data.frame(Node=unique(as.vector(outTab[,"ARG"])), Type="ARG")
nodeOut=rbind(lncNode, mrnaNode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)


ARGLncRNA=unique(as.vector(outTab[,"lncRNA"]))
ARGLncRNAexp=data[ARGLncRNA,]
ARGLncRNAexp=rbind(ID=colnames(ARGLncRNAexp), ARGLncRNAexp)
write.table(ARGLncRNAexp,file="ARGLncExp.txt",sep="\t",quote=F,col.names=F)
