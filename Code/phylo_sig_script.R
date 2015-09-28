library(phytools)
SB_tree <-read.newick("Raw/SB_Tree.phy")

SB_ids <- read.csv("Raw/SB_ids.csv", quote="\"")
SB_labels <- as.character(as.vector(SB_ids[,"fw_label"]))

mytree <- SB_tree
mytree$tip.label <- SB_labels

mytree <- collapse.singles(mytree)

plot(mytree)

datatree=treedata(tr,traits_out_df,sort=TRUE,warnings=TRUE)

max_rank = 50

Ps.out = rep("na",max_rank)
Ks.out = rep("na",max_rank)
Ks.out[1] = phylosig(datatree$phy,datatree$data[,1],method="K",test=TRUE)$K
Ps.out[1] = phylosig(datatree$phy,datatree$data[,1],method="K",test=TRUE)$P
for(i in 2:max_rank){
  temp = Test.Kmult(datatree$data[,1:i],datatree$phy)
  Ks.out[i] = temp$phy.signal[1,1]
  Ps.out[i] = temp$pvalue[1,1]
}

phy_test = datatree$phy
data_test = as.matrix(datatree$data)
data.r_test = as.matrix(data_test[sample(nrow(data_test)),])
rownames(data.r_test)<-rownames(data_test)

phylosig(phy_test,data_test,method="K",test=TRUE)
phylosig(phy_test,data.r_test,method="K",test=TRUE)



x<-as.matrix(datatree$data[,1:5])
N<-length(datatree$phy$tip.label)
ones<-array(1,N)
phy <- datatree$phy
C<-vcv.phylo(phy)
C<-C[row.names(x),row.names(x)]
a.obs<-colSums(solve(C))%*%x/sum(solve(C)) #evol.vcv code
distmat<-as.matrix(dist(rbind(as.matrix(x),a.obs)))
MSEobs.d<-sum(distmat[(1:N),(N+1)]^2) #sum distances root vs. tips
eigC <- eigen(C)
D.mat<-solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors))
dist.adj<-as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0)))
MSE.d<-sum(dist.adj[(1:N),(N+1)]^2) #sum distances for transformed data)
K.denom<-(sum(diag(C))- N*solve(t(ones)%*%solve(C)%*%ones)) / (N-1)
K.stat<-(MSEobs.d/MSE.d)/K.denom


par(mfcol=c(1,2))

Ks.out.spl <- smooth.spline(Ks.out)
plot(Ks.out)
lines(Ks.out.spl)

Ps.out.spl <- smooth.spline(Ps.out)
plot(Ps.out)
lines(Ps.out.spl)
abline(h=0.05,col="red")

Test.Kmult(datatree_in$data[,1:10],datatree_in$phy)

phylosig(datatree_in$phy,datatree_in$data[,1],method="K",test=TRUE)
contMap(mytree,traits_in[,1],outline=FALSE, fsize=0,legend=0, palette="heatmap",cex = 0.5)

<-fastBM(mytree,nsim=1000)
Y<-fastBM(mytree,nsim=1000)
K<-apply(X,2,Test.Kmult,phy=mytree,"lambda")
quantile(,c(0.05,0.95))


Test.Kmult(X[,c(1,2)],mytree)

for(i in 2:15){Test.Kmult(X[,1:i],mytree)}


par(fg="transparent")
plotTree(mytree,fsize=0.4,ylim=c(-1,length(mytree$tip.label)))
lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
ss<-sort(unique(gr))
colors<-setNames(cm.colors(14)[1:length(ss)],ss)
par(fg="black")
#add.simmap.legend(colors=colors,vertical=FALSE,x=0.25, y=-1,prompt=FALSE)
colors<-sapply(gr,function(gr,y) y[which(names(y)==gr)], y=colors)
tt<-gsub("_"," ",mytree$tip.label)
text(lastPP$xx[1:length(tt)],lastPP$yy[1:length(tt)], tt,cex=0.6,col=colors,pos=4,offset=0.1,font=2)

phenogram(mytree,scale(TraitOUT[,2]))


TraitsPCA_v <- TraitPCA$scores
names(TraitsPCA_v) <- names(TraitPCA$communality)


groups_ebd = c(13,13,13,13,8,13,13,13,13,13,11,13,13,13,13,9,13,13,13,13,13,13,9,13,13,9,13,13,9,13,13,13,14,13,13,13,11,14,13,13,13,13,13,13,13,13,13,12,11,9,11,13,11,11,13,9,13,13,9,11,13,9,8,13,13,13,13,13,13,13,13,13,13,13,9,11,13,8,11,9,9,11,8,8,9,13,9,8,13,9,9,8,8,11,11,7,9,9,9,11,13,9,7,9,11,11,11,9,11,9,11,9,10,13,13,9,13,11,11,11,13,7,8,9,11,11,13,13,11,3,3,6,4,4,3,4,4,4,3,3,6,3,6,3,3,6,3,3,4,2,2,1,1,1,2,2,1,1,5,5,6)
names(groups) <- mytree$tip


tree<-pbtree(n=25)
X<-fastBM(tree,nsim=2)
phylomorphospace(tree,X,xlab="trait 1",ylab="trait 2")
phylosig(mytree,gr,method="lambda",test=TRUE)
phylosig(mytree,gr,test=TRUE)


Trait <- Trait[order(rownames(Trait)),]
T_d <- dist(Trait, method = "euclidean")
B <- hclust(T_d, method="ward")
plot(B)
groups <- cutree(B, k=14)
rect.hclust(B, k=14, border="red")
cont <- table(groups,traits_cluster_13[,3])
matchClasses(cont,method="rowmax")
cont_out <- table(traits_cluster_13[,3],groups_out_W_M)
matchClasses(ccont,method="rowmax")

randIndex(cont_out)

clusters.kmeans <- kmeans(Trait, 14, algorithm="MacQueen")
plot(Trait, col=clusters.kmeans$cluster)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=clusters.kmeans$cluster, pch=clusters.kmeans$cluster)

T_Data <- as.data.frame(Trait)
T_PCA <- princomp(~.,data=T_Data,cor=TRUE)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=groups, pch=groups)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=traits_cluster_13[,3], pch=traits_cluster_13[,3])

