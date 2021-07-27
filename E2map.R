###这是画图之前的材料整理工作
setwd("/Users/xuebozhao/Documents/LuLab/Wheatpolyploidy/Sample/")
sampleall =  read.table("Data_File_S1.csv",head =T,sep = ",")
sample = read.table("AandAB_50.txt",head = F)
index = match(sample$V1,sampleall$ID)

sampleinfo = sampleall[index,]
write.table(sampleinfo, "sampleinfo.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


#######################总的map的图（总的分类信息和配色）
##打开总的文件
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/")
info1 = read.delim("map_file_r.csv",head =T,sep=",")
unicountry = unique(info1$Origin_country)
##打开分组文件
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group")
info = info1[!is.na(info1$Latitude),]
A = read.table("AA_taxa.txt",head =F)
AB = read.table("AABB_taxa.txt",head =F)
ABD = read.table("AABBDD_taxa.txt",head =F)
D = read.table("DD_taxa.txt",head =F)
S = read.table("SS_taxa.txt",head =F)
orderA = match(A[,1],info[,1])
orderAB = match(AB[,1],info[,1])
orderABD = match(ABD[,1],info[,1])
orderD = match(D[,1],info[,1])
orderS = match(S[,1],info[,1])
infoA = info[orderA,][!is.na(info[orderA,][,1]),]
infoAA = infoA[,(7:8)]
infoAB = info[orderAB,][!is.na(info[orderAB,][,1]),]
infoAABB = infoAB[,(7:8)]
infoABD = info[orderABD,][!is.na(info[orderABD,][,1]),]
infoAABBDD = infoABD[,(7:8)]
infoD = info[orderD,][!is.na(info[orderD,][,1]),]
infoDD = infoD[,(7:8)]
infoS = info[orderS,][!is.na(info[orderS,][,1]),]
infoSS = infoS[,(7:8)]
##开始啦
library(OpenStreetMap)
library(ggplot2)
map <- openmap(c(70,-179), c(-70,179),type='esri')
a = projectMercator(infoAA[,1],infoAA[,2])
ab = projectMercator(infoAABB[,1],infoAABB[,2])
abd = projectMercator(infoAABBDD[,1],infoAABBDD[,2])
d = projectMercator(infoDD[,1],infoDD[,2])
s = projectMercator(infoSS[,1],infoSS[,2])
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
pink = rgb(1,0.5,1,0.6)
skyblue = rgb(0,1,1,0.6)
green = rgb(152/255,251/255,152/255,0.8)
yellow = rgb(255/255,185/255,15/255,0.6)
plot(map,raster=TRUE)
points(a,col=red,pch = 19,cex = 2) 
points(ab,col=blue,pch = 15,cex = 2)
points(abd,col=pink,pch = 17,cex = 2)
points(d,col=yellow,pch = 19,cex = 2)
points(s,col=green,pch = 19,cex = 2)
legend(-18000000,150000,box.lty=0,bg="transparent",c("AA","SS","DD", "AABB","AABBDD"), 
       pch=c(19,19,19,15,17),col=c(red,green,yellow,blue,pink),cex=1.8,x.intersp=0.5,y.intersp=0.4)

###换一种画法
#install.packages("maps") ## 安装R包
#install.packages("ggplot2")
library("ggplot2") ### 加载ggplot2
library("maps") ### 加载地图包
world_map <- map_data("world") ### 导入地图数据 
head(world_map) ### 看看地图数据是啥样子
world_map$size <- 1:dim(world_map)[1] ### 生产一些随机的数据

### 画图开始 ####
ggplot(world_map, aes(x = long, y = lat, group = group, fill = size)) +
  geom_polygon(colour = "white",size=0.2) + scale_fill_distiller(palette = "YlGnBu",direction = 1)


###现在展示野生二粒的南北分群
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/")
info1 = read.delim("map_file_r.csv",head =T,sep=",")
info = info1[!is.na(info1$Latitude),]
unicountry = unique(info1$Origin_country)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression")
S = read.table("South_wildemmer.txt",head =F)
N = read.table("North_wildemmer.txt",head =F)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
Dom = read.table("sub_Domesticated_emmer.txt",head =F)
orderS = match(S[,1],info[,1])
orderN = match(N[,1],info[,1])
orderDom = match(Dom[,1],info[,1])
infoS = info[orderS,][!is.na(info[orderS,][,1]),]
infoSS = infoS[,(7:8)]
infoN = info[orderN,][!is.na(info[orderN,][,1]),]
infoNN = infoN[,(7:8)]
infoDom = info[orderDom,][!is.na(info[orderDom,][,1]),]
infoDomDom = infoDom[,(7:8)]
library(OpenStreetMap)
library(ggplot2)
dev.off()
map <- openmap(c(-2,-10), c(60,60),type='esri')
a = projectMercator(infoSS[,1],infoSS[,2])
ab = projectMercator(infoNN[,1],infoNN[,2])
dom = projectMercator(infoDomDom[,1],infoDomDom[,2])
plot(map,raster=TRUE)
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
skyblue = rgb(0,1,1,0.6)
yellow = rgb(168/255, 110/255, 87/255, 0.8) #a86e57
points(a,col=red,pch = 19,cex = 2) 
points(ab,col=blue,pch = 15,cex = 2)
points(dom,col=yellow,pch = 17,cex = 2)

##展示speltoids和wild emmer和dom emmer的地域分布情况，表示基因是发生在同域的地方
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/")
info1 = read.delim("map_file_r.csv",head =T,sep=",")
info = info1[!is.na(info1$Latitude),]
unicountry = unique(info1$Origin_country)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
S = read.table("sub_Speltoides.txt",head =F)
Wild = read.table("sub_Wild_emmer.txt",head =F)
Dom = read.table("sub_Domesticated_emmer.txt",head =F)
orderS = match(S[,1],info[,1])
orderWild = match(Wild[,1],Wild[,1])
orderDom = match(Dom[,1],info[,1])
infoS = info[orderS,][!is.na(info[orderS,][,1]),]
infoSS = infoS[,(7:8)]
infoWild = info[orderWild,][!is.na(info[orderWild,][,1]),]
infoWildWild = infoWild[,(7:8)]
infoDom = info[orderDom,][!is.na(info[orderDom,][,1]),]
infoDomDom = infoDom[,(7:8)]
library(OpenStreetMap)
library(ggplot2)
dev.off()
map <- openmap(c(-2,-10), c(60,60),type='esri')
a = projectMercator(infoSS[,1],infoSS[,2])
ab = projectMercator(infoWildWild[,1],infoWildWild[,2])
dom = projectMercator(infoDomDom[,1],infoDomDom[,2])
plot(map,raster=TRUE)
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
skyblue = rgb(0,1,1,0.6)
yellow = rgb(168/255, 110/255, 87/255, 0.8) #a86e57
points(a,col=red,pch = 19,cex = 2) 
points(ab,col=blue,pch = 15,cex = 2)
points(dom,col=yellow,pch = 17,cex = 2)
legend(2000000,3500000,box.lty=0,bg="transparent",c("Speltoids", "Wild emmer","Dometicated emmer"),
       pch=c(19,15,17),col=c(red,blue,yellow),cex=1.2,x.intersp=1,y.intersp=0.7)


#####PCA的结果展示
rm(list=ls())
dev.off()
###Alineage
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/Alineage")
a = read.table("A.eigenvec",head =T)
head(a[,1:4])
av = read.table("A.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/treeRangeColor_E2.txt",head =F)
dim(a)
dim(pca_color)
index = match(pca_color$V1,a$IID)
head(explained_ratio)
#par(bg="black")
#par(mfrow=c(2,1), oma=c(4.5,4, 4, 2.5), mar=c(0.5,1,0.5,1), cex=1)
#layout(matrix(c(1:2), 1, 2, byrow = F), widths=rep(1,2), heights=c(rep(1,2),0.05))
#bty="l" 只保留左和下两条边框
par(mfrow=c(1,2))
par(mfrow=c(2,1))
plot(a[index,3:4],xlab="PC1 (42.01%)",ylab="PC2 (13.70%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
#mtext("A lineage",side = 3 ,las = 0, line  = -1.4,cex = 1.5,at = -0.017)
plot(a[index,5:6],xlab="PC3 (3.73%)",ylab="PC4 (2.85%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
#mtext("A lineage",side = 3 ,las = 0, line  = 0.5,cex = 1.5,at = -0.24)

###Blineage
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/Blineage")
a = read.table("B.eigenvec",head =T)
head(a[,1:4])
av = read.table("B.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/treeRangeColor_E2.txt",head =F)
dim(a)
dim(pca_color)
index = match(pca_color$V1,a$IID)
head(explained_ratio)
par(mfrow=c(1,2))
plot(a[index,3:4],xlab="PC1 (32.51%)",ylab="PC2 (4.98%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
plot(a[index,5:6],xlab="PC3 (4.64%)",ylab="PC4 (3.95%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)


###Dlineage
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/Dlineage")
a = read.table("D.eigenvec",head =T)
head(a[,1:4])
av = read.table("D.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/treeRangeColor_E2.txt",head =F)
dim(a)
dim(pca_color)
index = match(pca_color$V1,a$IID)
head(explained_ratio)
par(mfrow=c(1,2))
plot(a[index,3:4],xlab="PC1 (55.38%)",ylab="PC2 (9.41%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
plot(a[index,5:6],xlab="PC3 (3.48%)",ylab="PC4 (3.13%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)


###Dlineage里面的ABD
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/Dlineage")
a = read.table("D.ABDtaxa.eigenvec",head =T)
head(a[,1:4])
av = read.table("D.ABDtaxa.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/treeRangeColor_E2_forABD.txt",head =F)
dim(a)
dim(pca_color)
index = match(pca_color$V1,a$IID)
head(explained_ratio)
#par(mfrow=c(2,1))
plot(a[index,3:4],xlab="PC1 (5.59%)",ylab="PC2 (3.88%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
plot(a[index,5:6],xlab="PC3 (3.57%)",ylab="PC4 (3。42%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)
ggplot(a[index,5:6],aes(x=PC3,y=PC4))+ 
  geom_point(aes(color=as.character(pca_color[,3])),size=5) + 
  theme_bw() +
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))
#######Alineage
rm(list=ls())
library(ape)
library(RColorBrewer)
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/sNMF/Alineage")
library(viridis)
col11=brewer.pal(12,'Set3')
headA = read.table("labelA.2.Q",head =F)
group = read.delim("A_order.txt",head=F)
index = match(group$V1,headA$V1)
K = 6
par(mfrow=c(K+1,1), oma=c(4.5,4, 4, 2.5), mar=c(0.3,1,0.3,1), cex=1)
layout(matrix(c(1:(K+1)), K+1, 1, byrow = TRUE), widths=rep(1,K), heights=c(rep(1,K),0.2))
for ( k in c(2:(K+1))){
  a= read.delim(paste("labelA.",k,".Q",sep=""),head =F)
  aa = a[,-1]
  barplot(t(as.matrix(aa[index,])),border=NA,width = 2,axisnames=F,col=col11,space = 0,axes = T,names.arg=NULL)
  mtext(c<-paste('K = ',k),side = 4 ,las = 1, line  = -1.7,cex = 1,at = 0.5)
}
mtext("Ancestry Fraction",side = 2 ,las = 0, line  = 2.6 ,cex = 1.2,at = 3.2)
a = 1
b = 486
plot(NULL, xlim=c(a-19, b+19), ylim=c(0, 20), main="",axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
polygon(c(a, a, b,b), c(20,3,3,20),col="white", border="grey")
rect(1,20,32,3,col="#FB8072",border = NA)
rect(32,20,63,3,col="#BEBADA",border = NA)
rect(63,20,92,3,col="#80B1D3",border = NA)

rect(92,20,142,3,col="#FDB462",border = NA)
rect(142,20,171,3,col="#B3DE69",border = NA)
rect(171,20,178,3,col="#8B636C",border = NA) 
rect(178,20,181,3,col="#FFB5C5",border = NA) 

rect(181,20,193,3,col="#FF1493",border = NA)
rect(193,20,203,3,col="#CD919E",border = NA)
rect(203,20,221,3,col="#FFFFB3",border = NA)
rect(221,20,231,3,col="#EEA9B8",border = NA)
rect(231,20,241,3,col="#FFE4E1",border = NA) 

rect(241,20,246,3,col="#9B30FF",border = NA) #Club_wheat
rect(241,20,251,3,col="#BC80BD",border = NA) #Macha
rect(251,20,275,3,col="#CDB38B",border = NA) #Spelt
rect(275,20,280,3,col="#FFE1FF",border = NA) #Indian_dwarf_wheat
rect(280,20,285,3,col="#87CEFF",border = NA) #Yunan_wheat
rect(285,20,290,3,col="#36648B",border = NA) #Tibetan_semi_wild
rect(290,20,295,3,col="#B452CD",border = NA) #Xinjiang_wheat
rect(295,20,300,3,col="#836FFF",border = NA) #Vavilovii

rect(300,20,426,3,col="#CCEBC5",border = NA) #Landrace
rect(426,20,485,3,col="#FFED6F",border = NA) #Cultivar
rect(485,20,487,3,col="#CDB38B",border = NA) #Synthetic


#######Blineage
rm(list=ls())
library(ape)
library(RColorBrewer)
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/sNMF/Blineage")
library(viridis)
col11=brewer.pal(12,'Set3')
headA = read.table("labelB.2.Q",head =F)
group = read.delim("B_order.txt",head=F)
index = match(group$V1,headA$V1)
K = 6
par(mfrow=c(K+1,1), oma=c(4.5,4, 4, 2.5), mar=c(0.3,1,0.3,1), cex=1)
layout(matrix(c(1:(K+1)), K+1, 1, byrow = TRUE), widths=rep(1,K), heights=c(rep(1,K),0.2))
for ( k in c(2:(K+1))){
  a= read.delim(paste("labelB.",k,".Q",sep=""),head =F)
  aa = a[,-1]
  barplot(t(as.matrix(aa[index,])),border=NA,width = 2,axisnames=F,col=col11,space = 0,axes = T,names.arg=NULL)
  mtext(c<-paste('K = ',k),side = 4 ,las = 1, line  = -1.7,cex = 1,at = 0.5)
}
mtext("Ancestry Fraction",side = 2 ,las = 0, line  = 2.6 ,cex = 1.2,at = 3.2)
a = 1
b = 405
plot(NULL, xlim=c(a-17, b+17), ylim=c(0, 20), main="",axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
polygon(c(a, a, b,b), c(20,3,3,20),col="white", border="grey")
rect(1,20,11,3,col="#779970",border = NA)

rect(92-81,20,142-81,3,col="#FDB462",border = NA)
rect(142-81,20,171-81,3,col="#B3DE69",border = NA)
rect(171-81,20,178-81,3,col="#8B636C",border = NA) 
rect(178-81,20,181-81,3,col="#FFB5C5",border = NA) 

rect(181-81,20,193-81,3,col="#FF1493",border = NA)
rect(193-81,20,203-81,3,col="#CD919E",border = NA)
rect(203-81,20,221-81,3,col="#FFFFB3",border = NA)
rect(221-81,20,231-81,3,col="#EEA9B8",border = NA)
rect(231-81,20,241-81,3,col="#FFE4E1",border = NA) 

rect(241-81,20,246-81,3,col="#9B30FF",border = NA) #Club_wheat
rect(241-81,20,251-81,3,col="#BC80BD",border = NA) #Macha
rect(251-81,20,275-81,3,col="#CDB38B",border = NA) #Spelt
rect(275-81,20,280-81,3,col="#FFE1FF",border = NA) #Indian_dwarf_wheat
rect(280-81,20,285-81,3,col="#87CEFF",border = NA) #Yunan_wheat
rect(285-81,20,290-81,3,col="#36648B",border = NA) #Tibetan_semi_wild
rect(290-81,20,295-81,3,col="#B452CD",border = NA) #Xinjiang_wheat
rect(295-81,20,300-81,3,col="#836FFF",border = NA) #Vavilovii

rect(300-81,20,426-81,3,col="#CCEBC5",border = NA) #Landrace
rect(426-81,20,485-81,3,col="#FFED6F",border = NA) #Cultivar
rect(485-81,20,487-81,3,col="#CDB38B",border = NA) #Synthetic


#######Dlineage
rm(list=ls())
library(ape)
library(RColorBrewer)
rm(list=ls())
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/sNMF/Dlineage")
library(viridis)
col11=brewer.pal(12,'Set3')
headA = read.table("labelD.2.Q",head =F)
group = read.delim("D_order.txt",head=F)
index = match(group$V1,headA$V1)
K = 6
par(mfrow=c(K+1,1), oma=c(4.5,4, 4, 2.5), mar=c(0.3,1,0.3,1), cex=1)
layout(matrix(c(1:(K+1)), K+1, 1, byrow = TRUE), widths=rep(1,K), heights=c(rep(1,K),0.2))
for ( k in c(2:(K+1))){
  a= read.delim(paste("labelD.",k,".Q",sep=""),head =F)
  aa = a[,-1]
  barplot(t(as.matrix(aa[index,])),border=NA,width = 2,axisnames=F,col=col11,space = 0,axes = T,names.arg=NULL)
  mtext(c<-paste('K = ',k),side = 4 ,las = 1, line  = -1.7,cex = 1,at = 0.5)
}
mtext("Ancestry Fraction",side = 2 ,las = 0, line  = 2.6 ,cex = 1.2,at = 3.2)
a = 1
b = 309
plot(NULL, xlim=c(a-13, b+13), ylim=c(0, 20), main="",axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
polygon(c(a, a, b,b), c(20,3,3,20),col="white", border="grey")
rect(1,20,29,3,col="#1F78B4",border = NA)
rect(29,20,64,3,col="#8DD3C7",border = NA)

rect(241-177,20,246-177,3,col="#9B30FF",border = NA) #Club_wheat
rect(241-177,20,251-177,3,col="#BC80BD",border = NA) #Macha
rect(251-177,20,275-177,3,col="#CDB38B",border = NA) #Spelt
rect(275-177,20,280-177,3,col="#FFE1FF",border = NA) #Indian_dwarf_wheat
rect(280-177,20,285-177,3,col="#87CEFF",border = NA) #Yunan_wheat
rect(285-177,20,290-177,3,col="#36648B",border = NA) #Tibetan_semi_wild
rect(290-177,20,295-177,3,col="#B452CD",border = NA) #Xinjiang_wheat
rect(295-177,20,300-177,3,col="#836FFF",border = NA) #Vavilovii

rect(300-177,20,426-177,3,col="#CCEBC5",border = NA) #Landrace
rect(426-177,20,485-177,3,col="#FFED6F",border = NA) #Cultivar
rect(485-177,20,487-177,3,col="#CDB38B",border = NA) #Synthetic


####现在把structure的结果plot到地图上面
library(tess3r)
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/sNMF/inner/Alineage")
A_SNMF = read.table("labelA.6.Q",header = F)
A_SNMF1 = A_SNMF[,2:7]
barplot(as.qmatrix(A_SNMF1), border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix")
info = read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",header = T)
index = match(A_SNMF$V1,info$ID)
lati = info[index,c(7,8,11)]
datafdmap = cbind(A_SNMF,lati)
datafdmap = na.omit(datafdmap)
#
my.colors <- c("lightblue","red","green","tomato", "blue","orange")
my.palette <- CreatePalette(my.colors, 9)
datamatrix = datafdmap[,2:7]
datacoord = cbind(datafdmap$Logititude,datafdmap$Latitude)
plot(as.qmatrix(datamatrix), datacoord, method = "map.max", interpol = FieldsKrigModel(200),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex = 0.6,col.palette = my.palette)
#points(datacoord,col = "red",pch = 19,cex = 0.9)
infoAABB =  datafdmap[(datafdmap$Genome == "AABB"),]
coord_4 = cbind(infoAABB$Logititude,infoAABB$Latitude)
infoAABBDD =  datafdmap[(datafdmap$Genome == "AABBDD"),]
coord_6 = cbind(infoAABBDD$Logititude,infoAABBDD$Latitude)
points(coord_6,col = "#8B3A3A",pch = 17,cex = 1.1)
points(coord_4,col = "#00FFFF",pch = 19,cex = 1)

###把不同的地理位置的的亚种画到基因组上面
library(tess3r)
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration")
A_SNMF = read.table("labelA.txt",header = F)
A_SNMF1 = A_SNMF[,2:3]
barplot(as.qmatrix(A_SNMF1), border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix")
info = read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",header = T)
index = match(A_SNMF$V1,info$ID)
lati = info[index,c(7,8,11)]
datafdmap = cbind(A_SNMF,lati)
datafdmap = na.omit(datafdmap)
#
my.colors <- c("lightblue","red","green","tomato", "blue","orange")
my.palette <- CreatePalette(my.colors, 9)
datamatrix = datafdmap[,2:3]
datacoord = cbind(datafdmap$Logititude,datafdmap$Latitude)
plot(as.qmatrix(datamatrix), datacoord, method = "map.max", interpol = FieldsKrigModel(200),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex = 0.6,col.palette = my.palette)
#points(datacoord,col = "red",pch = 19,cex = 0.9)
infoAABB =  datafdmap[(datafdmap$Genome == "AABB"),]
coord_4 = cbind(infoAABB$Logititude,infoAABB$Latitude)
infoAABBDD =  datafdmap[(datafdmap$Genome == "AABBDD"),]
coord_6 = cbind(infoAABBDD$Logititude,infoAABBDD$Latitude)
points(coord_6,col = "#8B3A3A",pch = 17,cex = 1.1)
points(coord_4,col = "#00FFFF",pch = 19,cex = 1)

####
library(fields)
library(RColorBrewer)
library(mapplots)
library (LEA)
library(raster)
library(dplyr)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
asc.raster='/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/paux/world.asc'
long.all = dt$lon
lat.all = dt$lat
pcoa1.all =  PCs[,1]
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)
coordinates = cbind.data.frame(long=long.all,lat=lat.all)
cell_value_min = NULL
cell_value_max = NULL
colpalette = colorRampPalette(c("red","cyan","orange","blue"))(100)
mapadd=T
pts.size = .2
pts.shape = 19
theta=10
cluster = data.frame(pcoa1.all)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(pcoa1.all),max(pcoa1.all)))
#414
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim")
dat414 = read.csv("accessions_414_bioclim.csv",head=T)
dat414$Logititude = as.numeric(as.character(dat414$Logititude))
data414 = dat414[,c(1,7,8,12:30)]
dt414 = na.omit(data414)
coord414 = cbind.data.frame(long=dt414$Logititude,lat=dt414$Latitude)
points(coord414,pch=pts.shape,cex=0.2) #12*12

####PCA的结果展示
library("FactoMineR")
library("factoextra")
rv_pca_data <- t(root_vegetative_data)
batch_pca <- root_vegetative_data[1,] %>% as.character()
dat.pca <- PCA(rv_pca_data[,-1], graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = batch_pca,
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE,
             legend.title = "Groups"
)

######这是展示MDS的结果
library(shape)
library(cluster)
library(geosphere)
library(RColorBrewer)
library(maps)
max_k = 12
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/AB/A")
taxaID = read.table("A_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
nams <- map("world", namesonly=TRUE, plot=FALSE)
#
gen_dist <- read.table("A_IN_matrix.mdist")
geo_dist <- as.matrix(distm(meta[,c(3,2)], fun=distHaversine))
gen_null <- sum(unlist(gen_dist))
geo_null <- sum(unlist(geo_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Paired"), max_k),nrow=max_k)
###### INIT VARIABLES
# get global MDS
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100
## Prep evaluation tools 
dgen_Sil <- c(NA)	# total silhouette in genetic space for total clustering
dgeo_Sil <- c(NA)	# total silhouette in geographic space for total clustering
sgen_Sil <- c(NA)	# total silhouette in genetic space for filtered clustering
sgeo_Sil <- c(NA)	# total silhouette in geographic space for filtered clustering
dgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by total clustering
dgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by total clustering
sgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by filtered clustering
sgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by filtered clustering
sind_Prop <- c(NA)	# proportion of individuals that are kept after silhoutte filtering
##### INTERATE THROUGH DIFFERENT K's 
k =5
for (k in 12:max_k)
{
  
  ### 1 ### CALCULATE EMPIRICAL SILHOUTTE SCORES FROM SIMULATIONS OF EVOLVING 50:50 POPULATIONS
  # simulate admixture
  clu <- pam(gen_dist, k, diss =T, keep.diss = F)
  
  # pre-filter data for bad fit indviduals
  
  list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))
  
  ## Iteratively estimate silhouette thresholds and filter out individuals
  old <- 0
  
  while ( length(list) != old )
  {	
    sil <- c()
    old <- length(list)
    flt_dist <- gen_dist[list,list]
    flu <- pam(flt_dist,  k, diss =T, keep.diss = F)
    
    # choose which two populations (medoids) are closest 
    meds <- as.matrix(flt_dist[flu$id.med,flu$id.med])
    diag(meds) <- 1
    
    sist <- apply(meds,1,which.min)	
    
    sil <- c()
    
    for (I in 1:k)
    {
      pop1 <- I
      pop2 <- sist[I] 
      
      
      ind1 <- which(flu$clustering == pop1)
      ind2 <- which(flu$clustering == pop2)
      ## calculate mean of empirical silhouette scores for simulated admixtures
      
      
      sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*1, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*1), 2, mean)
      sim_dist <- cbind(rbind(flt_dist, sim), c(sim,0))
      rownames(sim_dist)[nrow(sim_dist)] <- 5050
      
      sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*1.2, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*0.8), 2, mean)
      sim_dist <- cbind(rbind(sim_dist, c(sim, 0 )), c(sim,0,0))
      rownames(sim_dist)[nrow(sim_dist)] <- 6040
      
      sim <- apply(rbind(apply(flt_dist[which(flu$clustering == pop1),], 2, mean)*0.8, apply(flt_dist[which(flu$clustering == pop2),], 2, mean)*1.2), 2, mean)
      sim_dist <- cbind(rbind(sim_dist, c(sim, 0, 0)), c(sim, 0, 0, 0))
      rownames(sim_dist)[nrow(sim_dist)] <- 4060
      
      # get empirical Silhouette scores
      
      
      sdis1 <- apply(sim_dist[(nrow(sim_dist)-2):nrow(sim_dist), c(ind1,ind1)],1,mean)
      sdis2 <- apply(sim_dist[(nrow(sim_dist)-2):nrow(sim_dist), c(ind2,ind2)],1,mean)		
      
      silmax <- apply(rbind(sdis1,sdis2),2,max)
      
      sils <- abs(sdis1-sdis2)/silmax
      
      
      sil <- c(sil,max(sils))
      
    }
    print(paste("k: ",k,", Kept indviduals: ",old, sep =""))
    ### 2 ### RUN KMEDOID CLUSTERING AND FILTER IT FOR SILHOUETTE  SCORES
    
    # run clustering algorithm
    tlu <- pam(gen_dist[list,list], k, diss =T, keep.diss = F)
    
    list <- c()
    for (o in 1:k)
    {
      list <- c(list, as.numeric(row.names(tlu$silinfo$widths[which(tlu$silinfo$widths[,1] == o & tlu$silinfo$widths[,3] >= sil[o]),])) )
    }
    if (length(list) <= k)
    {
      break
    }
  }
  ### END of filtering iterations
  
  ### GET COLOUR PALETTES FOR THIS K
  
  
  pal <- as.character(col_pal[,k-1])
  dim <- paste(pal, 80, sep ="")
  
  # write tables for PLINK files filtering 
  
  pops <- c()
  for (K in 1:k)
  {
    poop <- names(which.max(table(meta[list[which(clu$clustering[list] == K)],4])))
    if (poop %in% pops)
    {
      poop <- paste(poop,"A",sep="")
    }
    pops <- c(pops, poop )
  }
  
  write.table(cbind(sort(as.character(meta[list,1])), sort(as.character(meta[list,1])) ),paste("list2keepK",k,".txt",sep=""),quote = F, sep = "\t",row.names = F, col.names = F)
  opps <- rep(NA, nrow(meta))
  opps[list] <- pops[clu$clustering[list]]
  write.table(cbind(opps,as.character(meta[,1]))[order(meta[,1]),], paste("list2annoK",k,".txt",sep=""),quote = F, sep = " ",row.names = F, col.names = F)
  
  
  pdf(paste("maps_",k,".pdf",sep=""))
  # plot MAP
  par(mfrow = c(1,1))
  map("world", regions = nams)
  points(meta[,3],meta[,2],col = rgb(0.8,0.8,0.8,0.5),pch =19,cex=0.6)
  points(meta[list,3],meta[list,2],col = dim[clu$clustering[list]],pch =19,cex=0.6)
  dev.off()
  
  pdf(paste("clust_",k,".pdf",sep=""))
  # plot MDS
  
  plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
  text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
  if (length(list) > 0 )
  {
    text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
  }
  #	points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)
  
  # plot PIES
  if (length(list) > 0 )
  {
    ftab <- table(as.data.frame(cbind(as.character(meta[,4]),clu$clustering)))
    stab <- table(as.data.frame(cbind(as.character(meta[list,4]),clu$clustering[list])))
    
    fsum <- apply(ftab,2,sum)
    ssum <- apply(stab,2,sum)
    
    fscal <- fsum/max(fsum)
    sscal <- ssum/max(fsum)
    
    allc <- c(brewer.pal(11,"RdYlBu"), brewer.pal(11,"PRGn"))
    pal2 <- allc[1:length(unique(meta[list,4]))]
    
    par(mfrow = c(3,2))
    for (i in 1:k)
    {
      to <- grepl(paste("^",i,"$",sep=""), colnames(stab))
      if (any(to))
      {	
        pie(stab[,which(to)], radius = sqrt(sscal[which(to)]))
        text(x=1,y=1,labels = paste(ssum[which(to)],"/",fsum[i],sep=""))
      } else {
        pie(c(0,1), col = NA, radius = 0.001)
        text(x=1,y=1,labels = paste("0/",fsum[i],sep=""))
      }
      plotcircle(sqrt(fscal[which(to)]))
    }
  }
  
  # sort and order individuals
  cent <- matrix(NA, k,2)
  for (o in 1:k)
  {
    aa = na.omit(meta[which(clu$clustering == o),2])
    cent[o,1] <- mean(aa)
    cent[o,2] <- mean(aa)
  }
  ind <- hclust(dist(cent))
  ord0 <- order(ind$order)[clu$silinfo$widths[,1]]
  ord1 <- as.numeric(substring(row.names(clu$silinfo$widths[order(ord0),]),2))
  
  # plot HEAT
  par(mfrow = c(1,1))
  image(as.matrix(gen_dist[ord1,ord1]),ylab = NA ,xlab = NA, xaxt = 'n', yaxt = 'n')
  par(mfrow = c(4,1))
  barplot(height = rep(10,dim(gen_dist)[1]), width = rep(1,dim(gen_dist)[1]), rep(1,dim(gen_dist)[1]), col = pal[clu$clustering[ord1]], border = NA, space = 0, yaxt = 'n', xaxt = 'n' )
  barplot(height = rep(10,dim(gen_dist)[1]), width = rep(1,dim(gen_dist)[1]), rep(1,dim(gen_dist)[1]), col = pal2[meta[ord1,4]], border = NA, space = 0, yaxt = 'n', xaxt = 'n')
  
  dev.off()
  
  ### Calculate evaluation statistics
  dgeo_Pil <- c()
  sgeo_Pil <- c()
  
  dgen_WSS <- c()
  dgeo_WSS <- c()
  sgen_WSS <- c()
  sgeo_WSS <- c()
  
  for (J in 1:k)
  {
    pop1 <- J
    pop2 <- sist[J]
    
    ind1 <- which(clu$clustering == pop1)
    ind2 <- which(clu$clustering == pop2)	
    
    for (S in ind1)
    {
      d1 <- mean( unlist(geo_dist[S,ind1])) 
      d2 <- mean( unlist(geo_dist[S,ind2]))
      dgeo_Pil <- c(dgeo_Pil, abs(d1-d2)/max(c(d1,d2)))
    }
    
    dgen_WSS <- c(dgen_WSS, sum(unlist(gen_dist[ind1,ind1])) )
    dgeo_WSS <- c(dgeo_WSS, sum(unlist(geo_dist[ind1,ind1])) )
    
    ind1 <- intersect(ind1,list)						# only filtered inds
    ind2 <- intersect(ind2,list)						# only filtered inds
    
    for (S in ind1)
    {
      d1 <- mean( unlist(geo_dist[S,ind1]))   
      d2 <- mean( unlist(geo_dist[S,ind2]))
      sgeo_Pil <- c(sgeo_Pil, abs(d1-d2)/max(c(d1,d2)))
    }
    
    sgen_WSS <- c(sgen_WSS, sum(unlist(gen_dist[ind1,ind1])) )
    sgeo_WSS <- c(sgeo_WSS, sum(unlist(geo_dist[ind1,ind1])) )		
  }
  
  dgen_Sil <- c(dgen_Sil, mean(clu$silinfo$widths[,3]))
  dgeo_Sil <- c(dgeo_Sil, mean(dgeo_Pil))
  
  dgen_TSS <- c(dgen_TSS, sum(dgen_WSS)/gen_null)
  dgeo_TSS <- c(dgeo_TSS, sum(dgeo_WSS)/geo_null)
  
  sgen_Sil <- c(sgen_Sil, mean(c(tlu$silinfo$widths[,3],rep(0,length(clu$clustering)-length(list))) ))
  sgeo_Sil <- c(sgeo_Sil, mean(c(sgeo_Pil,rep(0,length(clu$clustering)-length(list))), na.rm = T ))
  
  sgen_null <- sum(unlist(gen_dist[list,list]))
  sgeo_null <- sum(unlist(geo_dist[list,list]))
  
  sgen_TSS <- c(sgen_TSS, sum(sgen_WSS)/sgen_null)
  sgeo_TSS <- c(sgeo_TSS, sum(sgeo_WSS)/sgeo_null)
  
  
  sind_Prop <- c(sind_Prop,length(list)/nrow(meta))
  
}

pdf("WholeEval.pdf", useDingbats = F)
plot(dgen_Sil, type = 'l', col = rgb(0.1,0.8,0.3), ylim = c(0,1), xlab  = "K", ylab = "Fraction")
lines(dgeo_Sil, col = rgb(0.8,0.1,0.3))
lines(1-dgen_TSS, col = rgb(0.4,0.8,0.5), lty = 2)
lines(1-dgeo_TSS, col = rgb(0.8,0.5,0.4), lty = 2)
dev.off()

pdf("FilterEval.pdf", useDingbats = F)
plot(sgen_Sil, type = 'l', col = rgb(0.1,0.8,0.3), ylim = c(0,1), xlab  = "K", ylab = "Fraction")
lines(sgeo_Sil, col = rgb(0.8,0.1,0.3))
lines(1-sgen_TSS, col = rgb(0.4,0.8,0.5), lty = 2)
lines(1-sgeo_TSS, col = rgb(0.8,0.5,0.4), lty = 2)
lines(sind_Prop, col = rgb(0.2,0.2,0.8), lty = 3)
dev.off()

###最后发现K=7的时候是最合理的，开始重新画A的图
max_k = 7
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/AB/A")
taxaID = read.table("A_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
nams <- map("world", namesonly=TRUE, plot=FALSE)
gen_dist <- read.table("A_IN_matrix.mdist")
geo_dist <- as.matrix(distm(meta[,c(3,2)], fun=distHaversine))
gen_null <- sum(unlist(gen_dist))
geo_null <- sum(unlist(geo_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Paired"), max_k),nrow=max_k)
###### INIT VARIABLES
# get global MDS
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100
## Prep evaluation tools 
dgen_Sil <- c(NA)	# total silhouette in genetic space for total clustering
dgeo_Sil <- c(NA)	# total silhouette in geographic space for total clustering
sgen_Sil <- c(NA)	# total silhouette in genetic space for filtered clustering
sgeo_Sil <- c(NA)	# total silhouette in geographic space for filtered clustering
dgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by total clustering
dgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by total clustering
sgen_TSS <- c(NA)		# proportion of genetic sum of squares explained by filtered clustering
sgeo_TSS <- c(NA)		# proportion of geographic sum of squares explained by filtered clustering
sind_Prop <- c(NA)	# proportion of individuals that are kept after silhoutte filtering
##### INTERATE THROUGH DIFFERENT K's 
k =7
### 1 ### CALCULATE EMPIRICAL SILHOUTTE SCORES FROM SIMULATIONS OF EVOLVING 50:50 POPULATIONS
# simulate admixture
clu <- pam(gen_dist, k, diss =T, keep.diss = F) #if TRUE (default for dist or dissimilarity objects), then x will be considered as a dissimilarity matrix
# pre-filter data for bad fit indviduals
#list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))
list <- as.numeric(substring(row.names(clu$silinfo$widths),2))
####开始画图，画地图 plot MAP
par(mfrow = c(1,1))
map("world", regions = nams)
points(meta[,3],meta[,2],col = rgb(0.8,0.8,0.8,0.5),pch =19,cex=0.6)
points(meta[list,3],meta[list,2],col = dim[clu$clustering[list]],pch =19,cex=0.6)
#画MDS的分布图， plot MDS
dev.off()

plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
if (length(list) > 0 )
{
  text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
}
# points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)
# 画饼图，plot PIES
dev.off()
if (length(list) > 0 )
{
  ftab <- table(as.data.frame(cbind(as.character(meta[,4]),clu$clustering)))
  stab <- table(as.data.frame(cbind(as.character(meta[list,4]),clu$clustering[list])))
  
  fsum <- apply(ftab,2,sum)
  ssum <- apply(stab,2,sum)
  
  fscal <- fsum/max(fsum)
  sscal <- ssum/max(fsum)
  
  allc <- c(brewer.pal(11,"RdYlBu"), brewer.pal(11,"PRGn"))
  pal2 <- allc[1:length(unique(meta[list,4]))]
  
  par(mfrow = c(3,3))
  for (i in 1:k)
  {
    to <- grepl(paste("^",i,"$",sep=""), colnames(stab))
    if (any(to))
    {	
      pie(stab[,which(to)], radius = sqrt(sscal[which(to)]))
      text(x=1,y=1,labels = paste(ssum[which(to)],"/",fsum[i],sep=""))
    } else {
      pie(c(0,1), col = NA, radius = 0.001)
      text(x=1,y=1,labels = paste("0/",fsum[i],sep=""))
    }
    plotcircle(sqrt(fscal[which(to)]))
  }
}
###最后发现K=7的时候画B的图
pal <- as.character(col_pal[,k-1])
dim <- paste(pal, 80, sep ="")
max_k = 8
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/AB/B")
taxaID = read.table("B_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
nams <- map("world", namesonly=TRUE, plot=FALSE)
gen_dist <- read.table("B_IN_matrix.mdist")
geo_dist <- as.matrix(distm(meta[,c(3,2)], fun=distHaversine))
gen_null <- sum(unlist(gen_dist))
geo_null <- sum(unlist(geo_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Set2"), max_k),nrow=max_k)
pal <- as.character(col_pal[,k-1])
dim <- paste(pal, 80, sep ="")
###### INIT VARIABLES
# get global MDS
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100
k =7
### 1 ### CALCULATE EMPIRICAL SILHOUTTE SCORES FROM SIMULATIONS OF EVOLVING 50:50 POPULATIONS
clu <- pam(gen_dist, k, diss =T, keep.diss = F) #if TRUE (default for dist or dissimilarity objects), then x will be considered as a dissimilarity matrix
#list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))
list <- as.numeric(substring(row.names(clu$silinfo$widths),2))
####开始画图，画地图 plot MAP
par(mfrow = c(1,1))
map("world", regions = nams)
points(meta[,3],meta[,2],col = rgb(0.8,0.8,0.8,0.5),pch =19,cex=0.6)
points(meta[list,3],meta[list,2],col = dim[clu$clustering[list]],pch =19,cex=0.6)
#画MDS的分布图， plot MDS
dev.off()
plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
if (length(list) > 0 )
{
  text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
}
# points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)
# 画饼图，plot PIES
dev.off()
if (length(list) > 0 )
{
  ftab <- table(as.data.frame(cbind(as.character(meta[,4]),clu$clustering)))
  stab <- table(as.data.frame(cbind(as.character(meta[list,4]),clu$clustering[list])))
  
  fsum <- apply(ftab,2,sum)
  ssum <- apply(stab,2,sum)
  
  fscal <- fsum/max(fsum)
  sscal <- ssum/max(fsum)
  
  allc <- c(brewer.pal(11,"RdYlBu"), brewer.pal(11,"PRGn"))
  pal2 <- allc[1:length(unique(meta[list,4]))]
  
  par(mfrow = c(3,3))
  for (i in 1:k)
  {
    to <- grepl(paste("^",i,"$",sep=""), colnames(stab))
    if (any(to))
    {	
      pie(stab[,which(to)], radius = sqrt(sscal[which(to)]))
      text(x=1,y=1,labels = paste(ssum[which(to)],"/",fsum[i],sep=""))
    } else {
      pie(c(0,1), col = NA, radius = 0.001)
      text(x=1,y=1,labels = paste("0/",fsum[i],sep=""))
    }
    plotcircle(sqrt(fscal[which(to)]))
  }
}
###最后发现K=4的时候画D的图
pal <- as.character(col_pal[,k-1])
dim <- paste(pal, 80, sep ="")
max_k = 8
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/D")
taxaID = read.table("D_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
nams <- map("world", namesonly=TRUE, plot=FALSE)
gen_dist <- read.table("D_IN_matrix.mdist")
geo_dist <- as.matrix(distm(meta[,c(3,2)], fun=distHaversine))
gen_null <- sum(unlist(gen_dist))
geo_null <- sum(unlist(geo_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Set2"), max_k),nrow=max_k)
pal <- as.character(col_pal[,k-1])
dim <- paste(pal, 80, sep ="")
###### INIT VARIABLES
# get global MDS
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100
k =4
### 1 ### CALCULATE EMPIRICAL SILHOUTTE SCORES FROM SIMULATIONS OF EVOLVING 50:50 POPULATIONS
clu <- pam(gen_dist, k, diss =T, keep.diss = F) #if TRUE (default for dist or dissimilarity objects), then x will be considered as a dissimilarity matrix
#list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))
list <- as.numeric(substring(row.names(clu$silinfo$widths),2))
####开始画图，画地图 plot MAP
par(mfrow = c(1,1))
map("world", regions = nams)
points(meta[,3],meta[,2],col = rgb(0.8,0.8,0.8,0.5),pch =19,cex=0.6)
points(meta[list,3],meta[list,2],col = dim[clu$clustering[list]],pch =19,cex=0.6)
#画MDS的分布图， plot MDS
dev.off()
plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
if (length(list) > 0 )
{
  text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
}
# points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)
# 画饼图，plot PIES
dev.off()
if (length(list) > 0 )
{
  ftab <- table(as.data.frame(cbind(as.character(meta[,4]),clu$clustering)))
  stab <- table(as.data.frame(cbind(as.character(meta[list,4]),clu$clustering[list])))
  
  fsum <- apply(ftab,2,sum)
  ssum <- apply(stab,2,sum)
  
  fscal <- fsum/max(fsum)
  sscal <- ssum/max(fsum)
  
  allc <- c(brewer.pal(11,"RdYlBu"), brewer.pal(11,"PRGn"))
  pal2 <- allc[1:length(unique(meta[list,4]))]
  
  par(mfrow = c(3,3))
  for (i in 1:k)
  {
    to <- grepl(paste("^",i,"$",sep=""), colnames(stab))
    if (any(to))
    {	
      pie(stab[,which(to)], radius = sqrt(sscal[which(to)]))
      text(x=1,y=1,labels = paste(ssum[which(to)],"/",fsum[i],sep=""))
    } else {
      pie(c(0,1), col = NA, radius = 0.001)
      text(x=1,y=1,labels = paste("0/",fsum[i],sep=""))
    }
    plotcircle(sqrt(fscal[which(to)]))
  }
}


####现在要做的是画出三个地方的空的地图来
library(OpenStreetMap)
library(ggplot2)
##这是新月沃地附近的地图
map <- openmap(c(0,10), c(55,70),type='esri')
plot(map,raster=TRUE)
##这是欧洲附近的地图
map <- openmap(c(15,-12), c(70,45),type='esri')
plot(map,raster=TRUE)
##这是欧洲附近的地图
map <- openmap(c(10,35), c(55,150),type='esri')
plot(map,raster=TRUE)


#######现在做的是migration的图
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot","EEMS_D")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-axes-flipped"),longlat = FALSE,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-output-PNGs"),longlat = TRUE,plot.height = 8,plot.width = 7,res = 600,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-demes-and-edges"),longlat = TRUE,add.grid = TRUE,col.grid = "gray90",lwd.grid = 2,
           add.outline = TRUE,col.outline = "blue",lwd.outline = 5,add.demes = TRUE,col.demes = "red",pch.demes = 5,min.cex.demes = 0.5,max.cex.demes = 1.5,out.png = FALSE)
##
library(rgdal)
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/test_figout/out1",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-cartographic-projections"),
           longlat = TRUE,projection.in = projection_none,projection.out = projection_mercator,out.png = FALSE)
library("rworldmap")
library("rworldxtra")
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-geographic-map"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,add.map = TRUE,col.map = "black",lwd.map = 5,out.png = FALSE)
library("RColorBrewer")
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-new-eems-colors"),longlat = TRUE,
           eems.colors = brewer.pal(11, "RdBu"),out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_africa <- map_world[which(map_world@data$continent == "Africa"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, 
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Similarly we can add points, lines, labels, etc. ## Here is how to add a set of colored "labels" on top of the ## migration/diversity rates and the Africa map
coords <- matrix(c(-10, 10, 10, 10, 30, 0, 40, -10, 30, -20), ncol = 2, byrow = TRUE) 
colors <- c("red", "green", "blue", "purple", "orange") 
labels <- LETTERS[1:5]
coords_merc <- sp::spTransform(SpatialPoints(coords, CRS(projection_none)), CRS(projection_mercator))
## `coords_merc` is a SpatialPoints structure ## but we only need the coordinates themselves
coords_merc <- coords_merc@coords
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-labels-projected"),longlat = TRUE,projection.in = projection_none,projection.out = projection_mercator,
           m.plot.xy = { plot(map_africa, col = NA, add = TRUE);text(coords_merc, col = colors, pch = labels, font = 2); },
           q.plot.xy = { plot(map_africa, col = NA, add = TRUE);text(coords_merc, col = colors, pch = labels, font = 2); },out.png = FALSE)
dev.off()

#######现在做的是migration的图---A
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/A/A_ex_out")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/A/A_ex_out_plot","EEMS_A")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/A/A_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_africa <- map_world[which(map_world@data$continent == "Africa"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, 
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
dev.off()

#######现在做的是migration的图---B
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/B/B_ex_out")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/B/B_ex_out_plot","EEMS_B")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/B/B_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_africa <- map_world[which(map_world@data$continent == "Africa"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, 
           add.demes =TRUE,col.demes = "black",min.cex.demes = 0.5, max.cex.demes = 2,add.grid = TRUE,
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
dev.off()

#######现在做的是migration的图---D
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot","EEMS_D")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_Eurasia <- map_world[which(map_world@data$continent == "Eurasia"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, col.outline = "gray90",
           add.demes =TRUE,col.demes = "black",min.cex.demes = 0.5, max.cex.demes = 2,
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
dev.off()

#######现在做的是migration的图---all
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/all/all_ex_out")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/all/all_ex_out_plot","EEMS_all")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/B/B_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_africa <- map_world[which(map_world@data$continent == "Africa"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, 
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
dev.off()



##现在是展示Vmap1.1的数据
data <- data.frame(
  name=letters[1:10],
  value=c(30,31,30,395,10,395,11,10,40,246)
)

# Uniform color
barplot(height=data$value, names=data$name, 
        col=c(rep("#BA9BC9",4),rep("#1F97FF",2),rep("#FFCE6E",4)),border = '#ffffff',
        horiz=T, las=1,xlim = c(0,400)
)



############现在做的是A的
library(shape)
library(cluster)
library(geosphere)
library(RColorBrewer)
max_k = 7
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/AB/A")
taxaID = read.table("A_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
gen_dist <- read.table("A_IN_matrix.mdist")
gen_null <- sum(unlist(gen_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Paired"), max_k),nrow=max_k)
k = 7
clu <- pam(gen_dist, k, diss =T, keep.diss = F) #if TRUE (default for dist or dissimilarity objects), then x will be considered as a dissimilarity matrix
# pre-filter data for bad fit indviduals
#list <- as.numeric(substring(row.names(clu$silinfo$widths[which(clu$silinfo$widths[,3] > 0),]),2))
list <- as.numeric(substring(row.names(clu$silinfo$widths),2))
#画MDS的分布图， plot MDS
dev.off()

plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
#text(-glob$points[,2],-glob$points[,1], label = meta[,4], col = rgb(0.8,0.8,0.8,0.5), cex = 0.5)
if (length(list) > 0 )
{
  text(-glob$points[list,2],-glob$points[list,1], label = meta[list,4],col = pal[clu$clustering[list]], cex = 0.5)
}

##MDS的点的分类
plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)

points(glob$points[clu$id.med,1],glob$points[clu$id.med,2], col = pal, cex = 1, pch =8)










