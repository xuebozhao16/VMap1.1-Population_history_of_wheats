library(ggtree)
library(treeio)
library(ape)
library(tidyverse)
library(RColorBrewer)

######这是算的是整个小麦属的分歧时间15 Ma的时间
mu <- 1e-3
my.colors <- c("#000000","#FFCE6E","#338732","#D39117","#4BB7B7","#447499","#645084","#000000")
my.labels <- c("barley","A001","A003","A005","A012","A028","A066","A070","A071","A076","A085","B023","B024","B034","B037","B045","B114","B115","B116","B119","B122","FIN-L1","GEO-L1","IRN-L3","UZB-L1","HCM")
tree <- read.beast("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/beast/species1line1_4971.mcc2.tre")
ggtree(tree) + geom_text(aes(label=node))
tree <- groupOTU(tree, list(seq(2,6),seq(7,11),seq(12,16),seq(17,21),seq(22,26)))
p <- ggtree(tree, aes(color=group))
revts(p) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", color="blue", alpha=.2, size=3, branch.length="height") +
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=my.colors, labels=my.labels) +
  #scale_x_continuous(limits=c(-0.018, 0.0005), breaks=seq(0,3,0.5)*(-mu), minor_breaks=seq(0,3,0.25)*(-mu), labels=seq(0,3,0.5))
  scale_x_continuous(limits=c(-0.018, 0.0005),breaks=seq(0,15,2.5)*(-mu), minor_breaks=seq(0,15,2.5)*(-mu), labels=seq(0,15,2.5))

######这是算的是整个A lineage的分歧时间2.5 Ma的时间,这个时间就短很多了
tree <- read.beast("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/A/A_1line1_4971.mc2.tre")
ggtree(tree) + geom_text(aes(label=node))
tree <- groupOTU(tree, list(seq(2,6)))
p <- ggtree(tree, aes(color=group))
revts(p) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", color="blue", alpha=.2, size=3, branch.length="height") +
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=my.colors, labels=my.labels) +
  scale_x_continuous(limits=c(-0.018, 0.0005), breaks=seq(0,3,0.5)*(-mu), minor_breaks=seq(0,3,0.25)*(-mu), labels=seq(0,3,0.5))

######这是算的是整个B lineage的分歧时间2.5 Ma的时间,这个时间就短很多了
tree <- read.beast("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/B/B_1line1_4971.mc2.tre")
ggtree(tree) + geom_text(aes(label=node))
tree <- groupOTU(tree, list(seq(2,16)))
p <- ggtree(tree, aes(color=group))
revts(p) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", color="blue", alpha=.2, size=3, branch.length="height") +
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=my.colors, labels=my.labels) +
  scale_x_continuous(limits=c(-0.008, 0.0005), breaks=seq(0,3,0.5)*(-mu), minor_breaks=seq(0,3,0.25)*(-mu), labels=seq(0,3,0.5))

######这是算的是整个D lineage的分歧时间2.5 Ma的时间,这个时间就短很多了
tree <- read.beast("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/D/D_1line1_4971.mc2.tre")
ggtree(tree) + geom_text(aes(label=node))
tree <- groupOTU(tree, list(seq(2,12)))
p <- ggtree(tree, aes(color=group))
revts(p) +
  theme_tree2() +
  geom_range(range="height_0.95_HPD", color="blue", alpha=.2, size=3, branch.length="height") +
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=my.colors, labels=my.labels) +
  scale_x_continuous(limits=c(-0.008, 0.0005), breaks=seq(0,3,0.5)*(-mu), minor_breaks=seq(0,3,0.25)*(-mu), labels=seq(0,3,0.5))


######这是为了得到一致根的树的做法，使用ape包, tol = 1e-3,
##这是在服务器上面做的，204上面是可以安装上ape的
library(ape)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/beast")
DD = read.table("species1line1_4971.mcc.tre",header = F, sep = "\t",quote = "\"")
DDtree = DD[1,]
DD.phy <- read.tree(text = as.character(DDtree)) #defining a tree in bracket form
plot(DD.phy)
plot(chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
              node = "root", S = 1, tol = 1e-8,
              CV = FALSE, eval.max = 500, iter.max = 500))
x <-DD.phy$edge.length
for (i in x){
  y<-x*7999
}
DD.phy$edge.length <-y
l <- 10^(-1:6)
cv <- numeric(length(l))

for (i in 1:length(l)){
  cv[i] <- sum(attr(chronopl(DD.phy, lambda = l[i], CV=TRUE), "D2"))
}
plot(l, cv)
plot(chronoMPL(DD.phy))
###
allapetree <- matrix("NA", nrow = 100, ncol = 1)
#allapetree[1,1] = "aaa"
for(i in 1:100){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol =  1e-8,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  #allapetree[i,1]=res[2,]
  write.tree(res,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/test_mcmc/apeTree/",i,".tree",sep=""))
}





















