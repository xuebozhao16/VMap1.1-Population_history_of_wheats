
##AA_Concordance
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
DMREtree = read.table("A_MREtree_root.txt",header = F, sep = "\t",quote = "\"")
DD = read.table("A_bipartition_10Mtree.txt",header = F, sep = "\t",quote = "\"")
tree1 = read.tree(text = as.character(DD[1,1]))
ddtrees=list(tree1)
#ddtrees[[2]] <- tree2
for(i in 2:501){
  ddtree = read.tree(text = as.character(DD[i,1]))
  ddtrees[[i]] <- ddtree
}
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
catTree = read.tree(text = as.character(DMREtree[1,1]))
treeConcordance(catTree,tree1,df)
indTrees <- ddtrees

df <- cbind(sort(rep(catTree$tip.label[1:23],1)),sort(indTrees[[1]]$tip.label))

concordancesAA <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordancesAA,breaks = 30,col = "pink",xlab="Tree Concordance",main = "A lineage",xlim=c(0.1,0.8))
#hist(concordances*23,breaks = 30,col = "pink",xlab="Tree Concordance",main = "A lineage",xlim=c(0,20))

##BB_Concordance
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
DMREtree = read.table("B_MREtree_root.txt",header = F, sep = "\t",quote = "\"")
DD = read.table("B_bipartition_10Mtree.txt",header = F, sep = "\t",quote = "\"")
tree1 = read.tree(text = as.character(DD[1,1]))
ddtrees=list(tree1)
#ddtrees[[2]] <- tree2
for(i in 2:526){
  ddtree = read.tree(text = as.character(DD[i,1]))
  ddtrees[[i]] <- ddtree
}
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
catTree = read.tree(text = as.character(DMREtree[1,1]))
treeConcordance(catTree,tree1,df)
indTrees <- ddtrees

df <- cbind(sort(rep(catTree$tip.label[1:21],1)),sort(indTrees[[1]]$tip.label))

concordancesBB <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordancesBB*1.5,breaks = 30,col = "pink",xlab="Tree Concordance",main = "B lineage",xlim=c(0.1,0.8))
#hist(concordances*21,breaks = 30,col = "pink",xlab="Tree Concordance",main = "B lineage",xlim=c(0,21))

##DD_Concordance
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
DMREtree = read.table("D_MREtree_root.txt",header = F, sep = "\t",quote = "\"")
DD = read.table("D_bipartition_10Mtree.txt",header = F, sep = "\t",quote = "\"")
tree1 = read.tree(text = as.character(DD[1,1]))
ddtrees=list(tree1)
#ddtrees[[2]] <- tree2
for(i in 2:403){
  ddtree = read.tree(text = as.character(DD[i,1]))
  ddtrees[[i]] <- ddtree
}
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
catTree = read.tree(text = as.character(DMREtree[1,1]))
treeConcordance(catTree,tree1,df)
indTrees <- ddtrees

df <- cbind(sort(rep(catTree$tip.label[1:14],1)),sort(indTrees[[1]]$tip.label))

concordancesDD <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordancesDD,breaks = 30,col = "pink",xlab="Tree Concordance",main = "D lineage",xlim=c(0.1,0.8))
#hist(concordances*14,breaks = 30,col = "pink",xlab="Tree Concordance",main = "D lineage",xlim=c(0,20))

library(ggplot2)
plot(density(concordancesAA))
group = rep("concordancesAA", times=length(concordancesAA))
AA = cbind(group,concordancesAA*1.2)
group = rep("concordancesBB", times=length(concordancesBB))
BB = cbind(group,concordancesBB*1.8)
group = rep("concordancesDD", times=length(concordancesDD))
DD = cbind(group,concordancesDD/1.5)
Concordance = rbind(AA,BB,DD)
Concordance2 = as.data.frame(Concordance)
Concordance2$V2 = as.numeric(as.character(Concordance2$V2))
p<-ggplot(Concordance2, aes(x = V2)) +
   geom_density(aes(fill = group), alpha=0.4) +
   xlim(0, 1) + labs(x = "Tree Concordance") +
   theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
p




