###现在表示4倍体和6倍体在全球的分布情况和冬春性
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/map")
info1 = read.delim("GrowthHabit_4_6x_map.csv",head =T,sep=",")
#unicountry = unique(info1$Origin_country)
info2 = info1[!is.na(info1$Latitude),]
info = info2[!is.na(info2$GrowthHabit),]
##打开分组文件
AABB= info[(info$Ploidy=="AABB"),]
AB_Winter = AABB[(AABB$GrowthHabit=="Winter"),6:7]
AB_Spring = AABB[(AABB$GrowthHabit=="Spring"),6:7]
AB_Facultative = AABB[(AABB$GrowthHabit=="Facultative"),6:7]
#AABBDD= info[(info$Ploidy=="AABBDD"),]
AABBDD= info[(info$Common_name=="Landrace"),]
ABD_Winter = AABBDD[(AABBDD$GrowthHabit=="Winter"),6:7]
ABD_Spring = AABBDD[(AABBDD$GrowthHabit=="Spring"),6:7]
ABD_Facultative = AABBDD[(AABBDD$GrowthHabit=="Facultative"),6:7]
##开始啦
library(OpenStreetMap)
library(ggplot2)
map <- openmap(c(70,-179), c(-70,179),type='esri')
ab1 = projectMercator(AB_Winter[,1],AB_Winter[,2])
ab2 = projectMercator(AB_Spring[,1],AB_Spring[,2])
ab3 = projectMercator(AB_Facultative[,1],AB_Facultative[,2])
abd1 = projectMercator(ABD_Winter[,1],ABD_Winter[,2])
abd2 = projectMercator(ABD_Spring[,1],ABD_Spring[,2])
abd3 = projectMercator(ABD_Facultative[,1],ABD_Facultative[,2])
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
pink = rgb(1,0.5,1,0.6)
skyblue = rgb(0,1,1,0.6)
green = rgb(152/255,251/255,152/255,0.8)
yellow = rgb(255/255,185/255,15/255,0.6)
plot(map,raster=TRUE)
points(ab1,col=red,pch = 17,cex = 1.5) 
points(ab2,col=blue,pch = 17,cex = 1.5)
points(ab3,col=yellow,pch = 17,cex = 1.5)
points(abd1,col=red,pch = 19,cex = 1.5) 
points(abd2,col=blue,pch = 19,cex = 1.5)
points(abd3,col=yellow,pch = 19,cex = 1.5)
legend(-18000000,150000,box.lty=0,bg="transparent",c("AABB","AABBDD"), 
       pch=c(17,19),col="#000000",cex=1.5,x.intersp=0.5,y.intersp=0.4)
legend(-18000000,-2500000,box.lty=0,bg="transparent",c("Winter","Spring","Facultative"), 
       pch=c(15),col=c(red,blue,yellow),cex=1.5,x.intersp=0.5,y.intersp=0.4)



#表示OTUB1这个基因在全球的分布情况
#学习网址 http://blog.sciencenet.cn/blog-724080-793453.html
# install.packages("reshape")
# install.packages("rworldmap")
# install.packages("rworldxtra")
require(reshape)
require (rworldmap)
require(rworldxtra)
rice <- read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/haplogroup/riceLonLat.txt", head=T) #读入今天要用到的数据（请下载附件数据）
rice_reshape <- rice_reshape <- cast(rice,Latitude+Longitude~Sub.population) #对数据进行预处理
rice_reshape <- as.data.frame(rice_reshape)
mapPies(rice_reshape,nameX="Longitude",nameY="Latitude",nameZs=c('ADMIX','AROMATIC','IND','TEJ','TRJ'),mapRegion='world',symbolSize=1,barOrient='vert',oceanCol="#BBFFFF",landCol="#FFDAB9") #饼图绘制，根据自己的需求调整大小和颜色
#mapBars(rice_reshape,nameX="Longitude",nameY="Latitude",nameZs=c('ADMIX','AROMATIC','IND','TEJ','TRJ'),mapRegion='world',symbolSize=5,barOrient='vert',oceanCol="dodgerblue",landCol="lightgreen") #柱状图绘制，根据自己的需求调整大小和颜色

require(reshape)
require (rworldmap)
require(rworldxtra)
wheat1 <- read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/haplogroup/wheatLonLat_norn.txt", head=T,sep="\t") 
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,nameX="Longitude",nameY="Latitude",nameZs=c('Winter','Spring','Facultative'),mapRegion='world',symbolSize=1,
        zColours=c(red,blue,yellow),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9") 
                                           
##表示这个基因的PCA的数据分析
red =  rgb(1,0,0,0.6)
blue = rgb(0,0,1,0.6)
yellow = rgb(255/255,185/255,15/255,0.6)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/PCA")
a = read.table("A.eigenvec",head =T)
head(a[,1:4])
max(a$PC1)
min(a$PC1)
max(a$PC2)
min(a$PC2)
max(a$PC3)
min(a$PC3)
max(a$PC4)
min(a$PC4)
av = read.table("A.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
head(explained_ratio)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/map")
info1 = read.delim("GrowthHabit_4_6x_map.csv",head =T,sep=",")
#unicountry = unique(info1$Origin_country)
info2 = info1[!is.na(info1$Latitude),]
info = info2[!is.na(info2$GrowthHabit),]
##打开分组文件
AABB= info[(info$Ploidy=="AABB"),]
AB_Winter = AABB[(AABB$GrowthHabit=="Winter"),1]
AB_Spring = AABB[(AABB$GrowthHabit=="Spring"),1]
AB_Facultative = AABB[(AABB$GrowthHabit=="Facultative"),1]
#AABBDD= info[(info$Ploidy=="AABBDD"),]
AABBDD= info[(info$Common_name=="Landrace"),]
ABD_Winter = AABBDD[(AABBDD$GrowthHabit=="Winter"),1]
ABD_Spring = AABBDD[(AABBDD$GrowthHabit=="Spring"),1]
ABD_Facultative = AABBDD[(AABBDD$GrowthHabit=="Facultative"),1]
index = match(AB_Winter,a$IID)
par(mfrow=c(1,2))
plot(a[match(AB_Winter,a$IID),3:4],xlab="PC1 (30.79%)",ylab="PC2 (20.18%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,
     col= red ,pch=17,bg="transparent",cex=1.5,ylim=c(-0.1,0.1),xlim=c(-0.1,0.1))
points(a[match(AB_Spring,a$IID),3:4],cex=1.5,col= blue ,pch=17)
points(a[match(AB_Facultative,a$IID),3:4],cex=1.5,col= yellow ,pch=17)
points(a[match(ABD_Winter,a$IID),3:4],cex=1.5,col= red ,pch=19)
points(a[match(ABD_Spring,a$IID),3:4],cex=1.5,col= blue ,pch=19)
points(a[match(ABD_Facultative,a$IID),3:4],cex=1.5,col= yellow ,pch=19)
plot(a[match(AB_Winter,a$IID),5:6],xlab="PC3 (12.11%)",ylab="PC4 (11.11%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,
     col= red ,pch=17,bg="transparent",cex=1.5,ylim=c(-0.03,0.015),xlim=c(-0.025,0.06))
points(a[match(AB_Spring,a$IID),5:6],cex=1.5,col= blue ,pch=17)
points(a[match(AB_Facultative,a$IID),5:6],cex=1.5,col= yellow ,pch=17)
points(a[match(ABD_Winter,a$IID),5:6],cex=1.5,col= red ,pch=19)
points(a[match(ABD_Spring,a$IID),5:6],cex=1.5,col= blue ,pch=19)
points(a[match(ABD_Facultative,a$IID),5:6],cex=1.5,col= yellow ,pch=19)
##Blineage
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/PCA")
a = read.table("B.eigenvec",head =T)
head(a[,1:4])
max(a$PC1)
min(a$PC1)
max(a$PC2)
min(a$PC2)
max(a$PC3)
min(a$PC3)
max(a$PC4)
min(a$PC4)
av = read.table("B.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
head(explained_ratio)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/map")
info1 = read.delim("GrowthHabit_4_6x_map.csv",head =T,sep=",")
#unicountry = unique(info1$Origin_country)
info2 = info1[!is.na(info1$Latitude),]
info = info2[!is.na(info2$GrowthHabit),]
##打开分组文件
AABB= info[(info$Ploidy=="AABB"),]
AB_Winter = AABB[(AABB$GrowthHabit=="Winter"),1]
AB_Spring = AABB[(AABB$GrowthHabit=="Spring"),1]
AB_Facultative = AABB[(AABB$GrowthHabit=="Facultative"),1]
#AABBDD= info[(info$Ploidy=="AABBDD"),]
AABBDD= info[(info$Common_name=="Landrace"),]
ABD_Winter = AABBDD[(AABBDD$GrowthHabit=="Winter"),1]
ABD_Spring = AABBDD[(AABBDD$GrowthHabit=="Spring"),1]
ABD_Facultative = AABBDD[(AABBDD$GrowthHabit=="Facultative"),1]
index = match(AB_Winter,a$IID)
par(mfrow=c(1,2))
plot(a[match(AB_Winter,a$IID),3:4],xlab="PC1 (32.59%)",ylab="PC2 (16.90%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,
     col= red ,pch=17,bg="transparent",cex=1.5,ylim=c(-0.02,0.06),xlim=c(-0.025,0.01))
points(a[match(AB_Spring,a$IID),3:4],cex=1.5,col= blue ,pch=17)
points(a[match(AB_Facultative,a$IID),3:4],cex=1.5,col= yellow ,pch=17)
points(a[match(ABD_Winter,a$IID),3:4],cex=1.5,col= red ,pch=19)
points(a[match(ABD_Spring,a$IID),3:4],cex=1.5,col= blue ,pch=19)
points(a[match(ABD_Facultative,a$IID),3:4],cex=1.5,col= yellow ,pch=19)
plot(a[match(AB_Winter,a$IID),5:6],xlab="PC3 (11.64%)",ylab="PC4 (10.09%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,
     col= red ,pch=17,bg="transparent",cex=1.5,ylim=c(-0.2,0.15),xlim=c(-0.2,0.1))
points(a[match(AB_Spring,a$IID),5:6],cex=1.5,col= blue ,pch=17)
points(a[match(AB_Facultative,a$IID),5:6],cex=1.5,col= yellow ,pch=17)
points(a[match(ABD_Winter,a$IID),5:6],cex=1.5,col= red ,pch=19)
points(a[match(ABD_Spring,a$IID),5:6],cex=1.5,col= blue ,pch=19)
points(a[match(ABD_Facultative,a$IID),5:6],cex=1.5,col= yellow ,pch=19)

###heatmap聚类
###OTUB1单倍型图谱
#教程https://www.jianshu.com/p/d86e4afe1065
#install.packages("pheatmap")
cols = c("#FFF5EB","#80B1D3","#FCCDE5")
library(heatmap.plus)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1/haplogroup")
##A
chr26Q_heatmap = read.table("chr37_TraesCS7A02G263900_2K_4_6.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
#annotation_c=info1[,c(1,10)]
annotation_c=info1[,c(10,12)]
rownames(annotation_c) <- colnames(data)
pheatmap(data,cluster_rows = F,cutree_cols = 2,color = cols,fontsize=4.5,annotation_col = annotation_c)
##B
chr26Q_heatmap = read.table("chr39_TraesCS7B02G161900_2K_4_6.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
#annotation_c=info1[,c(1,10)]
annotation_c=info1[,c(10,12)]
rownames(annotation_c) <- colnames(data)
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,annotation_col = annotation_c)
##D
chr26Q_heatmap = read.table("chr41_TraesCS7D02G264800_2K_6.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
#annotation_c=info1[,c(1,10)]
annotation_c=info1[,c(10,12)]
rownames(annotation_c) <- colnames(data)
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,annotation_col = annotation_c)


###和傅老师合作的，312个驯化相关的overlap基因里面，XPCLR值最高的一个gene，TraesCS5A02G233700的单倍型图谱
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/SNP")
chr26Q_heatmap = read.table("chr25_TraesCS5A02G233700_110K.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
heatmap(t(data),na.rm = T,scale = "none",Colv = NA,Rowv = NA, cexRow=0.5)

##这个是这个gene上游2K，下游2K的snp
cols = c("#FFF5EB","#80B1D3","#FCCDE5")
chr26Q_heatmap = read.table("chr25_TraesCS5A02G233700.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
heatmap(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.5)
#install.packages("heatmap.plus")
library(heatmap.plus)
heatmap.plus(t(data),col=cols,scale = "none",Colv = NA,Rowv = NA, cexRow=0.5)

##ggseqlogo 绘制序列分析图
#https://www.plob.org/article/12630.html
#install.packages("ggseqlogo")
library(ggplot2)
library(ggseqlogo)
data(ggseqlogo_sample)
fasta = "ACGTATG
ATGTATG
ACGTATG
ACATATG
ACGTACG"
ggplot()+geom_logo(seqs_dna$MA0004.1)+theme_logo()
ggplot()+geom_logo(fasta)+theme_logo()


##单倍型聚类图
cols = c("#FFF5EB","#80B1D3","#FCCDE5")
library(heatmap.plus)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/SNP")
##A
chr26Q_heatmap = read.table("chr25_TraesCS5A02G233700.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
info1 = read.delim("lodging_2_4x_map.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(10,12)]
rownames(annotation_c) <- colnames(data)
pheatmap(data,cluster_rows = F,cutree_cols = 6,color = cols,fontsize=4.5,annotation_col = annotation_c)

##查看TraesCS5A02G233700这个基因在驯化这一个过程中的变化，即野生一粒和栽培一粒，野生二粒和栽培二粒，去掉倍性之间的影响
##野生一粒和栽培一粒
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/SNP")
chr26Q_heatmap = read.table("mac_WD_A.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
info1 = read.delim("lodging_2_4x_map.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(10,12)]
rownames(annotation_c) <- colnames(data)
labels_row = chr26Q_heatmap$POS
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,
         annotation_col = annotation_c,labels_row = labels_row)  + geom_seqlogo(data)
library(ggplot2)
library(ggseqlogo)
fasta_file = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/SNP/mac_WD_A.fa")
fasta = fasta_file[,2]
data(ggseqlogo_sample)
ggplot()+geom_logo(seqs_dna$MA0004.1)+theme_logo()
#CS	CATGAGATCAGTCACATGGGTGTGAAGTGATAGCAG
ggplot()+geom_logo(as.character(fasta),method = "probability",stack_width = 0.7)+theme_logo()
ggplot() +geom_logo(as.character("CATGAGATCAGTCACATGGGTGTGAAGTGATAGCAG"),method = "probability",stack_width = 0.7)+theme_logo()
##野生二粒和栽培二粒
chr26Q_heatmap = read.table("mac_WD_AB.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
# heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
# heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
info1 = read.delim("lodging_2_4x_map.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(10,12)]
rownames(annotation_c) <- colnames(data)
labels_row = chr26Q_heatmap$POS
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,
         annotation_col = annotation_c,labels_row = labels_row)  + geom_seqlogo(data)
fasta_file = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/SNP/mac_WD_AB.fa")
fasta = fasta_file[,2]
ggplot()+geom_logo(as.character(fasta),method = "probability",stack_width = 0.7)+theme_logo()
#CS	CTGGGCGTCGCTCAGAAGTCGT
ggplot()+geom_logo(as.character("CTGGGCGTCGCTCAGAAGTCGT"))+theme_logo()

###表达量
##TraesCS5A02G233700 TraesCS5B02G232200 TraesCS5D02G240600
rm(list=ls())
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/RNA")
par(mfrow=c(3,1),oma=c(4.5,4, 4, 2.5), mar=c(2,1,1,1), cex=1)
layout(matrix(c(1:3), 3, 1, byrow = TRUE), widths=rep(1,3), heights=c(rep(1,3),0.002))
gene7A=read.delim("expression_7A.csv",head=T,sep = ",")
data1 = gene7A[,2]
data2= gene7A[,3]
aa=barplot(as.matrix(t(data1)),ylab = "7A expression Level(TPM)",col = "lavender",ylim = c(0,4))
arrows(aa,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7A expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 2)
gene7B=read.delim("expression_7B.csv",head=T,sep = ",")
data1 = gene7B[,2]
data2= gene7B[,3]
bb=barplot(as.matrix(t(data1)),ylab = "7B expression Level(TPM)",col = "DarkSeaGreen1",ylim = c(0,4))
arrows(bb,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7B expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 2)
gene7D=read.delim("expression_7D.csv",head=T,sep = ",")
data1 = gene7D[,2]
data2= gene7D[,3]
dd=barplot(as.matrix(t(data1)),ylab = "7D expression Level(TPM)",col = "LightGoldenrod3",ylim = c(0,4))
arrows(dd,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7D expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 2)
cap1 = gene7D[,1]
text(dd, y=par()$usr[3]-0.2,srt=30, xpd=T, adj=1, labels = cap1,cex = 1)

#####查看TraesCS5A02G233700这个基因在驯化这一个过程中的变化,换了一个新的包
library(ggmsa)
library(ggplot2)
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
ggmsa(protein_sequences, 164, 213, color = "Chemistry_AA")
nt_sequence <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
ggmsa(nt_sequence, color = "Chemistry_NT")
##结论，不好使

###画这个基因的各个部分的分类图
a = 449254000
b = 449260000
plot(NULL, xlim=c(a, b), ylim=c(0, 5), main="",axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
polygon(c(a, a, b,b), c(2.7,3,3,2.7),col="grey", border="grey")
rect(449254955,2.7,449255062,3,col="#00868B",border = NA)
rect(449255063,2.7,449255198,3,col="#00FFFF",border = NA)
rect(449255307,2.7,449255709,3,col="#00FFFF",border = NA)
rect(449256551,2.7,449257487,3,col="#00FFFF",border = NA)
rect(449259163,2.7,449259241,3,col="#00868B",border = NA)

##查看TraesCS5A02G233700这个基因在驯化这一个过程中的变化，即野生一粒和栽培一粒，野生二粒和栽培二粒，去掉倍性之间的影响
##野生一粒和栽培一粒indel
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/Crispr_gene/TraesCS5A02G233700/indel")
chr26Q_heatmap = read.table("mac_WD_Aindel.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
info1 = read.delim("lodging_2_4x_map.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(10,12)]
rownames(annotation_c) <- colnames(data)
labels_row = chr26Q_heatmap$POS
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,
         annotation_col = annotation_c,labels_row = labels_row)  + geom_seqlogo(data)
##野生二粒和栽培二粒inel
chr26Q_heatmap = read.table("mac_WD_ABindel.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
#heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
# heatmap.plus(t(data),scale = "none",col=cols,Colv = NA)
# heatmap.plus((data),scale = "none",col=cols,Rowv = NA)
library(pheatmap)
info1 = read.delim("lodging_2_4x_map.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(10,12)]
rownames(annotation_c) <- colnames(data)
labels_row = chr26Q_heatmap$POS
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,
         annotation_col = annotation_c,labels_row = labels_row)  + geom_seqlogo(data)
