#install.packages("igraph")


library(igraph)               
nodefile="net.node.txt"       
edgefile="net.network.txt"    
outfile="network.pdf"         
lncRNAcol="#00AFBB"           
NRGcol="#FC4E07"              
setwd("C:\\6.network")     


node.data=read.table(nodefile, header=T, sep="\t", check.names=F)
edge.data=read.table(edgefile, header=T, sep="\t", check.names=F)
color=ifelse(node.data$Type=="lncRNA", lncRNAcol, NRGcol)
value=ifelse(node.data$Type=="lncRNA", 2, 5)
fontSize=ifelse(node.data$Type=="lncRNA", 0.01, 0.65)
node=data.frame(id=node.data$Node,label=node.data$Node,color=color,shape="dot",value=value,fontSize=fontSize)
edge=data.frame(from=edge.data$lncRNA,to=edge.data$NRG,length=100,arrows="middle",smooth=TRUE,shadow=FALSE,weight=edge.data$cor)


d=data.frame(p1=edge$from, p2=edge$to, weight=abs(edge$weight))
g=graph.data.frame(d,directed = FALSE)
E(g)$color="grey"
V(g)$size=node$value[match(names(components(g)$membership),node$label)]
V(g)$shape="sphere"
V(g)$lable.cex=node$fontSize[match(names(components(g)$membership),node$label)]
V(g)$color=node$color[match(names(components(g)$membership),node$label)]


pdf(outfile, width=8, height=7)
layout(mat=matrix(c(1,2,1,2),nc=2), height=c(1,11))
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
legend('center',legend=c('lncRNA','NRG'),col=c(lncRNAcol,NRGcol),pch=16,bty="n",ncol=2,cex=2)
vertex.frame.color = node$color
edge_col=E(g)$color
plot(g,layout=layout_on_sphere,vertex.size=V(g)$size,vertex.label=node$label,vertex.label.cex=V(g)$lable.cex,edge.width =0.05,edge.arrow.size=0,vertex.label.color=NULL,vertex.frame.color=NA,edge.color=edge_col,vertex.label.font=2)
dev.off()
