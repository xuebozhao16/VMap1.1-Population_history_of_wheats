library(ggplot2)
library(ggcorrplot)
library(ggthemes)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/majorAllele/mac3_10M/")
data = read.table("Dtree_10M_topologyNUM.txt",head = F)
data1 = (data[,-1])
data2 = t(data1)
corr = round(cor(data2),3)  ##这是计算数据两两之间的相关性
p.mat = round(cor_pmat(data1),6)  #检验库之间的相关性的P值
ggcorrplot(corr)
ggcorrplot(corr,method = "circle")
ggcorrplot(corr,method = "circle",hc.order = TRUE,hc.method = "ward.D",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#6D9EC1","white","#db0d0d"),lab_size = 0.05)  ##这是进行聚类
data = read.table("Btree_10M_topologyNUM.txt",head = F)
data1 = (data[,-1])
data2 = t(data1)
corr = round(cor(data2),3)  ##这是计算数据两两之间的相关性
p.mat = round(cor_pmat(data1),6)  #检验库之间的相关性的P值
ggcorrplot(corr)
ggcorrplot(corr,method = "circle")
ggcorrplot(corr,method = "circle",hc.order = TRUE,hc.method = "ward.D",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#6D9EC1","white","#db0d0d"),lab_size = 0.05)  ##这是进行聚类





require(ape)
phy <- rtree(n=20)
plot(phy)
plot(phy,type="fan")
chronopl(phy, lambda=0, age.min = 1, age.max = NULL,
         node = "root", S = 1, tol = 1e-8,
         CV = FALSE, eval.max = 500, iter.max = 500)
plot(chronopl(phy, lambda=0, age.min = 1, age.max = NULL,
              node = "root", S = 1, tol = 1e-8,
              CV = FALSE, eval.max = 500, iter.max = 500))

DD.phy <- read.tree(text = "((Meyeri:0.00134539227411050941,(((Landrace:0.00224065264812887965,(((Cultivar:0.02269980609340177188,Club_wheat:0.00612447058578657715):0.00789993087666918808,Macha:0.00450723668837937734):0.00166840907104719616,((Indian_dwarf_wheat:0.00000100000050002909,(Tibetan_semi_wild:0.00000100000050002909,(Yunan_wheat:0.00000100000050002909,Xinjiang_wheat:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909,(Spelt:0.00759640474201685786,Vavilovii:0.00000100000050002909):0.00243645562193337863):0.00000100000050002909):0.00028715979690966575):0.00274292475451662224,Synthetic:0.00437696662506303642):0.02049431969156784855,Strangulata:0.00000100000050002909):0.00744946847804108317):0.00036736431507447353,Anathera:0.00036736431507447353);
") #defining a tree in bracket form
plot(DD.phy)
plot(chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
              node = "root", S = 1, tol = 1e-8,
              CV = FALSE, eval.max = 500, iter.max = 500))
##
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/assembly/D")
DD = read.table("D_best_30.txt",header = F, sep = "\t",quote = "\"")
DDtree = DD[1,]
DD.phy <- read.tree(text = as.character(DDtree)) #defining a tree in bracket form
plot(DD.phy)
plot(chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
              node = "root", S = 1, tol = 1e-8,
              CV = FALSE, eval.max = 500, iter.max = 500))
###
allapetree <- matrix("NA", nrow = 100, ncol = 1)
#allapetree[1,1] = "aaa"
for(i in 1:100){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol = 1e-8,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  #allapetree[i,1]=res[2,]
  write.tree(res,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/assembly/D/apeTree/",i,".tree",sep=""))
}

##现在是做10M的A B D的
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
#DD = read.table("A_best_10Mtree.txt",header = F, sep = "\t",quote = "\"")
#DD = read.table("B_best_10Mtree.txt",header = F, sep = "\t",quote = "\"")
DD = read.table("D_best_10Mtree.txt",header = F, sep = "\t",quote = "\"")
for(i in 1:403){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol = 1e-8,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  write.tree(res,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/",i,".tree",sep=""))
}
 

###这是10M的window里面的树的比较
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("D_MREtree",header = F, sep = "\t",quote = "\"")
AA = read.table("D_best_10Mtreechr.txt",header = F, sep = "\t",quote = "\"")
a <- read.tree(text = as.character(AMREtree[1,]))
b <- read.tree(text = as.character(AA[1,4]))
comparePhylo(a, b, plot = T, force.rooted = TRUE)
res = comparePhylo(a, b, plot = F, force.rooted = TRUE)
res1 = res$messages[8]
res11 = strsplit(res1,split = " ")
#substr(res1,1,2)
res111 = as.numeric(res11[[1]][1])
res2 = paste(AA[1,1],AA[1,2],AA[1,3],res111,sep="\t")
##DD
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("D_MREtree",header = F, sep = "\t",quote = "\"")
AA = read.table("D_best_10Mtreechr.txt",header = F, sep = "\t",quote = "\"")
for(i in 1:403){
  bb = read.tree(text = as.character(AA[i,4]))
  res = comparePhylo(a, bb, plot = F, force.rooted = TRUE)
  res1 = res$messages[8]
  res11 = strsplit(res1,split = " ")
  res111 = as.numeric(res11[[1]][1])
  res2 = paste(AA[i,1],AA[i,2],AA[i,3],res111,sep="\t")
  write(res2,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score/",i,".compare.txt",sep=""))
}
##AA
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("A_MREtree",header = F, sep = "\t",quote = "\"")
AA = read.table("A_best_10Mtreechr.txt",header = F, sep = "\t",quote = "\"")
for(i in 1:501){
  bb = read.tree(text = as.character(AA[i,4]))
  res = comparePhylo(a, bb, plot = F, force.rooted = TRUE)
  res1 = res$messages[8]
  res11 = strsplit(res1,split = " ")
  res111 = as.numeric(res11[[1]][1])
  res2 = paste(AA[i,1],AA[i,2],AA[i,3],res111,sep="\t")
  write(res2,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score/",i,".compare.txt",sep=""))
}
##BB
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("B_MREtree",header = F, sep = "\t",quote = "\"")
AA = read.table("B_best_10Mtreechr.txt",header = F, sep = "\t",quote = "\"")
for(i in 1:526){
  bb = read.tree(text = as.character(AA[i,4]))
  res = comparePhylo(a, bb, plot = F, force.rooted = TRUE)
  res1 = res$messages[8]
  res11 = strsplit(res1,split = " ")
  res111 = as.numeric(res11[[1]][1])
  res2 = paste(AA[i,1],AA[i,2],AA[i,3],res111,sep="\t")
  write(res2,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score/",i,".compare.txt",sep=""))
}

###画图啦
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score")
##AA
DDscore = read.table("A_treeScore_sort.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="A lineage 10M window",xlab="Tree Score",xlim = c(0,23))
polygon(density(DDscore$V4), col="#9BCD9B",border="#9BCD9B")
plot(density(DDscore$V4/23),main="A lineage 10M window",xlab="Tree Score/clades",xlim = c(0,1))
polygon(density(DDscore$V4/23), col="#FFE4B5",border="#FFE4B5")
aa = hist(DDscore$V4)
aa$counts
##BB
DDscore = read.table("B_treeScore_sort.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="B lineage 10M window",xlab="Tree Score",xlim = c(0,21))
polygon(density(DDscore$V4), col="#9BCD9B",border="#9BCD9B")
plot(density(DDscore$V4/23),main="B lineage 10M window",xlab="Tree Score/clades",xlim = c(0,1))
polygon(density(DDscore$V4/23), col="#FFE4B5",border="#FFE4B5")
aa = hist(DDscore$V4)
aa$counts
##DD
DDscore = read.table("D_treeScore_sort.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="D lineage 10M window",xlab="Tree Score",xlim = c(0,14))
polygon(density(DDscore$V4), col="#9BCD9B",border="#9BCD9B")
plot(density(DDscore$V4/14),main="D lineage 10M window",xlab="Tree Score/clades",xlim = c(0,1))
polygon(density(DDscore$V4/14), col="#FFE4B5",border="#FFE4B5")
aa = hist(DDscore$V4)
aa$counts

##D的图
library(RColorBrewer)
colDD = brewer.pal(6,'Set3')
source('~/Documents/LuLab/wheatSpeciation/geneTree/DensityPlot_wheat.R', echo=TRUE)
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/centromerer.txt",head=T)
##DD
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score")
S = read.table("D_treeScore_chr.txt",header = F, sep = "\t",quote = "\"")
head(S)
max(S$V3)
min(S$V3)
Densitplot_wheat(S,cent,legend.min = 6,legend.max = 12,legend.len = 12,col = c("darkgreen", "yellow", "red"),legend.x.intersp=1)
##AA
S = read.table("A_treeScore_chr.txt",header = F, sep = "\t",quote = "\"")
head(S)
max(S$V3)
min(S$V3)
Densitplot_wheat(S,cent,legend.min = 6,legend.max = 12,legend.len = 12,col = c("darkgreen", "yellow", "red"),legend.x.intersp=1)
##BB
S = read.table("B_treeScore_chr.txt",header = F, sep = "\t",quote = "\"")
head(S)
max(S$V3)
min(S$V3)
Densitplot_wheat(S,cent,legend.min = 6,legend.max = 12,legend.len = 12,col = c("darkgreen", "yellow", "red"),legend.x.intersp=1)


###test
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("D_MREtree",header = F, sep = "\t",quote = "\"")
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score/test")
AA = read.table("test3.txt",header = F, sep = "\t",quote = "\"")
for(i in 1:95){
  bb = read.tree(text = as.character(AA[i,4]))
  res = comparePhylo(a, bb, plot = F, force.rooted = TRUE)
  res1 = res$messages[8]
  res11 = strsplit(res1,split = " ")
  res111 = as.numeric(res11[[1]][1])
  res2 = paste(AA[i,1],AA[i,2],AA[i,3],res111,sep="\t")
  write(res2,paste("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree/tree_score/test/",i,".compare.txt",sep=""))
}
par(mfrow=c(3,1), oma=c(1,2, 2, 2), mar=c(3,1,1,1), cex=1)
DDscore = read.table("test1_score.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="",xlab="Tree Score",xlim = c(0,14))
polygon(density(DDscore$V4), col="#FB8072",border="#FB8072")
DDscore = read.table("test2_score.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="",xlab="Tree Score",xlim = c(0,14))
polygon(density(DDscore$V4), col="#80B1D3",border="#80B1D3")
DDscore = read.table("test3_score.txt",header = F, sep = "\t",quote = "\"")
plot(density(DDscore$V4),main="",xlab="Tree Score",xlim = c(0,14))
polygon(density(DDscore$V4), col="#FDB462",border="#FDB462")

######introgression 6A
dev.off()
cols = c("#FFF5EB","#FCCDE5","#fa75b9")
library(heatmap.plus)
library(pheatmap)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/6A")
chr26Q_heatmap = read.table("260M_430M_200_sort.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/introgression/Dstatistic/haplo_2A5B")
info1 = read.delim("accession_list_forR.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(2,3)]
rownames(annotation_c) <- colnames(data)
ann_colors = list(
  Ploidy = c(AABB = "#b3afb1", AABBDD = "#615f60"),
  Common_name = c(wild_emmer = "#527f9e", domesticated_emmer = "#9c86bf", free_threshing_tetraploid = "#00F5FF",
                  landrace_EA="#FDD662",landrace_EU="#8F9E76",landrace_WA="#FF8C8C")
)
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,annotation_col = annotation_c,
         annotation_colors = ann_colors)

####Q gene的单倍型的分布
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
require(reshape)
require (rworldmap)
require(rworldxtra)
wheat1 <- read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/haploNet/wheatLonLat_nornQ.txt", head=T,sep="\t") 
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,nameX="Longitude",nameY="Latitude",nameZs=c('q','Q'),mapRegion='world',symbolSize=1,
        zColours=c(red,blue),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9") 


####继续测试树的比较方法
#https://cran.r-project.org/web/packages/treespace/vignettes/tipCategories.html 这是学习的代码
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
AMREtree = read.table("D_MREtree",header = F, sep = "\t",quote = "\"")
AA = read.table("D_best_10Mtreechr.txt",header = F, sep = "\t",quote = "\"")
a <- read.tree(text = as.character(AMREtree[1,]))
b <- read.tree(text = as.character(AA[1,4]))
comparePhylo(a, b, plot = T, force.rooted = TRUE)
all.equal(a, b)
identical(a, b)
# load treespace and packages for plotting:
library(treespace)
library(RColorBrewer) 
library(ggplot2)
library(reshape2)
# set colour scheme
pal <- brewer.pal(3,"Dark2")

suppressWarnings(RNGversion("3.5.0"))
set.seed(948)
# set colour scheme
pal2 <- brewer.pal(8,"Dark2")

# create a "base" (category-level) tree
baseTree <- rtree(8)
baseTree$tip.label <- letters[8:1]

tree1  <- read.tree(text = as.character(AA[1,4]))
tree2 <- read.tree(text = as.character(AA[2,4]))
tree3 <- read.tree(text = as.character(AA[3,4]))
tree4 <- read.tree(text = as.character(AA[4,4]))
tree5 <- read.tree(text = as.character(AA[5,4]))
tree6 <- read.tree(text = as.character(AA[6,4]))
tipsMRCAdepths(tree1)

# set up colour palettes
tipcolors1 <- c(rep(pal2[[1]],3),rep(pal2[[2]],3),rep(pal2[[3]],3),rep(pal2[[4]],3),rep(pal2[[5]],3),rep(pal2[[6]],3),rep(pal2[[7]],3),rep(pal2[[8]],3))
tipcolors4 <- c(rep(pal2[[1]],4),rep(pal2[[2]],4),rep(pal2[[3]],4),rep(pal2[[4]],4),rep(pal2[[5]],4),rep(pal2[[6]],4),rep(pal2[[7]],4),rep(pal2[[8]],4)) # colours for 4 tips
tipcolors6 <- c(rep(pal2[[1]],6),rep(pal2[[2]],6),rep(pal2[[3]],6),rep(pal2[[4]],6),rep(pal2[[5]],6),rep(pal2[[6]],6),rep(pal2[[7]],6),rep(pal2[[8]],6)) # colours for 6 tips
# prepare tip colours for plotting
tree1TipOrder <- sapply(tree1$tip.label, function(x) which(tree1$tip.label==x))
layout(matrix(c(1,4,2,5,3,6), 2,3))
plot(tree1, tip.color=tipcolors1, no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=2)
plot(tree2, tip.color=tipcolors1, no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=2)
plot(tree3, tip.color=tipcolors1, no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.8)
plot(tree4, tip.color=tipcolors1[tree1TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)
plot(tree5, tip.color=tipcolors1[tree1TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)
plot(tree6, tip.color=tipcolors1[tree1TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)
##
trees <- list(tree1,tree2,tree3,tree4,tree5,tree6)
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
dists <- relatedTreeDist(trees,df)
dists
dists <- as.matrix(dists)
colnames(dists) <- rownames(dists) <- c("Tree 1", "Tree 2", "Tree 3", "Tree 4",
                                        "Tree 5", "Tree 6")
melted_dists <- melt(dists, na.rm=TRUE)

ggheatmap <- ggplot(data = melted_dists, aes(Var2, Var1, fill = value))+
  geom_tile(color = "darkgrey")+
  scale_fill_gradient2(low = "white", high = "firebrick2",
                       name="Tree distance") +
  theme_minimal() + coord_fixed()

ggheatmap +
  geom_text(aes(Var2, Var1, label = signif(value,2)), color = "black", size = 8) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),
    legend.position = "none" )

###这是和一个树去做比较
catTree <- read.tree(text = as.character(AA[1,4]))
indTree1 <- read.tree(text = as.character(AA[2,4]))
indTree2 <- read.tree(text = as.character(AA[3,4]))
indTree3 <- read.tree(text = as.character(AA[4,4]))

plot(catTree, tip.color=pal,
     edge.width = 4, type="cladogram",
     label.offset= 0.5, font=4,
     edge.color=c(pal[[1]],"black",pal[[2]],pal[[3]]))
##
dev.off()
layout(matrix(1:3,1,3))
plot(indTree1, tip.color=c(rep(pal[[1]],4),rep(pal[[2]],2),rep(pal[[3]],3)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[1]],6),"black",rep(pal[[2]],3),rep(pal[[3]],5)))
plot(indTree2, tip.color=c(rep(pal[[1]],4),pal[[2]],rep(pal[[3]],2),pal[[2]],pal[[3]]),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[1]],6),rep("black",2),pal[[2]],pal[[3]],
                  rep("black",2),pal[[3]],pal[[2]],pal[[3]]))
plot(indTree3, tip.color=c(rep(pal[[3]],3),pal[[2]],rep(pal[[1]],2),pal[[2]],rep(pal[[1]],2)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[3]],4),rep("black",2),pal[[2]],pal[[1]],
                  rep("black",2),pal[[1]],pal[[2]],rep(pal[[1]],3)))
###
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
treeConcordance(catTree,indTree1,df)
treeConcordance(catTree,indTree2,df)
treeConcordance(catTree,indTree3,df)
n <- 5
reps <- 10
reftree <- rtree(n, tip.label=catTree$tip.label[1:n])
indTrees <- lapply(rep(seq(0,100,20),reps), function(x)
  simulateIndTree(reftree,itips=n,permuteTips=TRUE,tipPercent=x))

df <- cbind(sort(rep(catTree$tip.label[1:n],n)),sort(indTrees[[1]]$tip.label))

concordances <- sapply(indTrees, function(x) treeConcordance(reftree,x,df))

resultsTab <- as.data.frame(cbind(rep(seq(0,100,20),reps),concordances))
colnames(resultsTab) <- c("Percentage","Concordance")
resultsTab[,1] <- factor(resultsTab[,1], levels=seq(0,100,20))
plot <- ggplot(resultsTab, aes(x=Percentage, y=Concordance))

plot + geom_boxplot(aes(colour=Percentage)) + theme_bw() + guides(colour=FALSE) +
  xlab("Percentage of tips permuted") + ylim(c(0,1)) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18))

######开始用自己的数据去跑啦
library(treespace)
library(RColorBrewer) 
library(ggplot2)
library(reshape2)
##DD
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/genomebed/10M_tree")
DD = read.table("D_bipartition_10Mtree.txt",header = F, sep = "\t",quote = "\"")
tree1 = read.tree(text = as.character(DD[1,1]))
ddtrees=list(tree1)
#ddtrees[[2]] <- tree2
for(i in 2:403){
  ddtree = read.tree(text = as.character(DD[i,1]))
  ddtrees[[i]] <- ddtree
}
#trees <- list(tree1,tree2,tree3,tree4,tree5,tree6)
df <- cbind(sort(tree1$tip.label),sort(tree1$tip.label))
dists <- relatedTreeDist(ddtrees,df)
dists
dists <- as.matrix(dists)
melted_dists <- melt(dists, na.rm=TRUE)

ggheatmap <- ggplot(data = melted_dists, aes(Var2, Var1, fill = value))+
  geom_tile(color = "darkgrey")+
  scale_fill_gradient2(low = "white", high = "firebrick2",
                       name="Tree distance") +
  theme_minimal() + coord_fixed()
ggheatmap

ggheatmap +
  geom_text(aes(Var2, Var1, label = signif(value,2)), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),
    legend.position = "none" )

layout(matrix(1:3,1,3))
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

concordances <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordances,breaks = 30,col = "pink",xlab="Tree Concordance",main = "D lineage",xlim=c(0.1,0.8))
#hist(concordances*14,breaks = 30,col = "pink",xlab="Tree Concordance",main = "D lineage",xlim=c(0,20))

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

concordances <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordances,breaks = 30,col = "pink",xlab="Tree Concordance",main = "A lineage",xlim=c(0.1,0.8))
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

concordances <- sapply(indTrees, function(x) treeConcordance(catTree,x,df))
hist(concordances*1.5,breaks = 30,col = "pink",xlab="Tree Concordance",main = "B lineage",xlim=c(0.1,0.8))
#hist(concordances*21,breaks = 30,col = "pink",xlab="Tree Concordance",main = "B lineage",xlim=c(0,21))


########
#devtools::install_github("bbanbury/phrynomics")
library(phrynomics)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Cytoscape")
snpsold <- ReadSNP("mac_Bexontree.phy")
snps1 <- RemoveInvariantSites(snpsold)
snps2 <- RemoveNonBinary(snps1)
WriteSNP(snps2, file="BInvariantSite.phy", format="phylip")


#####计算cytoscape所需要的distance and edge
dev.off()
rm(list=ls())
#http://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Phylogenetic-trees.nb.html#network_to_cytoscape
library(ape)
library(phytools)
library(igraph)
library(RCy3)
library(ggplot2)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Cytoscape/B")
tree1 = read.table("RAxML_bestTree.Bexontree_1.tree",header = F, sep = "\t",quote = "\"")
tree <- phytools::read.newick(text=as.character(tree1[1,1]))
#tree <- phytools::read.newick(system.file("extdata","phylotree.newick", package="RCy3"))
ig0 <- ape::as.igraph.phylo(tree, FALSE) # boolean for whether tree is rooted or not
ig=mst(ig0, weights = NULL, algorithm = NULL)
tree2 = read.table("RAxML_bestTree.Bexontree_2.tree",header = F, sep = "\t",quote = "\"")
tree2 <- phytools::read.newick(text=as.character(tree2[1,1]))
ig1 <- ape::as.igraph.phylo(tree2, FALSE) 
ig11=mst(ig1, weights = NULL, algorithm = NULL)
#igraph::consensus_tree(ig, hrg = NULL, start = FALSE, num.samples = 2)
aa = cluster_edge_betweenness(ig, weights = E(ig)$weight,
                         directed = TRUE, edge.betweenness = TRUE, merges = TRUE,
                         bridges = TRUE, modularity = TRUE, membership = TRUE)
ig <- set_edge_attr(ig,'distance', value=tree$edge.length) # set distances as edge attributes
createNetworkFromIgraph(ig, title="phylotree", collection = "phylotree")
layoutNetwork(paste('force-directed',
                    'defaultEdgeWeight=3',
                    'defaultSpringCoefficient=5E-5',
                    'defaultSpringLength=80',
                    sep = ' '))
createColumnFilter('junctions', 'id', "^Node\\\\d+$", "REGEX")
junctions<-getSelectedNodes()
setNodeWidthBypass(junctions,1)
setNodeHeightBypass(junctions,1)
setNodeLabelBypass(junctions, "")
setEdgeLabelMapping('distance')
layoutNetwork(paste('force-directed',
                    'edgeAttribute="distance"',
                    'type="1 - normalized value"',
                    'defaultSpringCoefficient=5E-4',
                    'defaultSpringLength=50',
                    sep = ' '))

#####
require(stats)
X <- matrix(runif(200), 20, 10)
d <- dist(X)
PC <- prcomp(X)
M <- ape::mst(d)
opar <- par()
par(mfcol = c(2, 2))
plot(M)
plot(M, graph = "nsca")
plot(M, x1 = PC$x[, 1], x2 = PC$x[, 2])

###
library("ggnetwork")
load("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Cytoscape/mst/data/dist2009c.RData")
country09 = attr(dist2009c, "Label")
mstree2009 = ape::mst(dist2009c)
gr09 = graph.adjacency(mstree2009, mode = "undirected")
#gg = ggnetwork(gr09, arrow.gap = 0, layout = "fruchtermanreingold")
gg = ggnetwork(gr09,arrow.gap = 0)
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +  theme_blank() +
  #geom_nodetext(aes(label = name), color = "black", size = 2.5) +
  theme(plot.margin = unit(c(0, 1, 1, 6), "cm"))+
  guides(color = guide_legend(keyheight = 0.09, keywidth = 0.09,
                              title = "Countries")) + theme(legend.position = c(0, 0.14),
                                                            legend.background = element_blank(),
                                                            legend.text = element_text(size = 7))
# gg = ggnetwork(ig)
# ggplot(ig, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
#   geom_nodes(aes(color = vertex.attributes()), size = 2) +  theme_blank() +
#   theme(plot.margin = unit(c(0, 1, 1, 6), "cm"))+
#   guides(color = guide_legend(keyheight = 0.09, keywidth = 0.09,
#                               title = "Countries")) + theme(legend.position = c(0, 0.14),
#                                                             legend.background = element_blank(),
#                                                             legend.text = element_text(size = 7))

library("rworldmap")
mat = match(country09, countriesLow$NAME)
coords2009 = data.frame(
  lat = countriesLow$LAT[mat],
  lon = countriesLow$LON[mat],
  country = country09)
layoutCoordinates = cbind(
  x = jitter(coords2009$lon, amount = 15),
  y = jitter(coords2009$lat, amount = 8))
labc = names(table(country09)[which(table(country09) > 1)])
matc = match(labc, countriesLow$NAME)
dfc = data.frame(
  latc = countriesLow$LAT[matc],
  lonc = countriesLow$LON[matc],
  labc)
dfctrans = dfc
dfctrans[, 1] = (dfc[,1] + 31) / 93
dfctrans[, 2] = (dfc[,2] + 105) / 238
ggeo09 = ggnetwork(gr09, arrow.gap = 0, layout = layoutCoordinates)
#ggeo09 = ggnetwork(ig)
ggplot(ggeo09, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +
  theme_blank() +
  geom_label(data = dfctrans, aes(x = lonc, xend = lonc, y = latc, yend = latc,
                                  label = labc, fill = labc), colour = "white", alpha = 0.5, size = 3) +
  theme(legend.position = "none")


##
library(ggnetwork)
library(network)
library(sna)
n <- network(rgraph(10, tprob = 0.2), directed = FALSE)
ggdf=ggnetwork(n, layout = "fruchtermanreingold", cell.jitter = 0.75)
datf = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Cytoscape/mst/data/string_graph.txt", header = TRUE)
grs = graph_from_data_frame(datf[, c("node1", "node2")], directed = FALSE)
E(grs)$weight = 1
V(grs)$size = centralization.degree(grs)$res
#ggdf = ggnetwork(grs)
ggdf = fortify(grs, layout = layout.spring(grs), cell.jitter = 0)
ggplot(ggdf, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", curvature = 0.1, size = 0.95, alpha = 0.8)+
  geom_nodes(aes(x = x, y = y), size = 3, alpha = 0.5, color = "orange") +
  geom_nodelabel_repel(aes(label = name), size=4, color="#8856a7") +
  #  geom_nodetext(aes(label = vertex.names), size = 4, color = "#8856a7") +
  theme_blank() + theme(legend.position = "none")
par(oldpar)

cap1 = c("Speltoides","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat",
         "Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat",
         "Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace")
###
#install.packages("phyloseq")
require(phyloseq)
ps1  = readRDS("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Cytoscape/mst/data/ps1.rds")
sampledata = data.frame( sample_data(ps1))
d1 = as.matrix(phyloseq::distance(ps1, method="jaccard"))
gr = graph.adjacency(d1,  mode = "undirected", weighted = TRUE)
net = igraph::mst(gr)
V(net)$id = sampledata[names(V(net)), "host_subject_id"]
V(net)$litter = sampledata[names(V(net)), "family_relationship"]
gnet=ggnetwork(net)
ggplot(gnet, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter)) + theme_blank()+
  theme(legend.position="bottom")
       