######introgression 6A
dev.off()
cols = c("#FFF5EB","#FCCDE5","#fa75b9")
library(heatmap.plus)
library(pheatmap)
######这是测试的第一个
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block/")
chr26Q_heatmap = read.table("Haplo_chr5_61_77M_shuf1000_sort.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block")
info1 = read.delim("accession_list_ABD_forR.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(2,3)]
rownames(annotation_c) <- colnames(data)
ann_colors = list(
  Class = c(DD = "#615f60", ABD_orphon = "#a19d9f",Bread_wheat = "#d9d4d7"),
  Common_name = c(Anathera = "#527f9e", Meyeri = "#9c86bf", Strangulata = "#00F5FF",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F")
)
cols = c("#FFF5EB","#FCCDE5","#fa75b9")
pheatmap(data,cluster_rows = F,cutree_cols = 7,color = cols,fontsize=4.5,annotation_col = annotation_c,
         annotation_colors = ann_colors)

######这是测试的第二个
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block/")
chr26Q_heatmap = read.table("Haplo_chr5_104_105M_ABD.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block")
info1 = read.delim("accession_list_ABD_forR.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(2,3)]
rownames(annotation_c) <- colnames(data)
ann_colors = list(
  Class = c(DD = "#615f60", ABD_orphon = "#a19d9f",Bread_wheat = "#d9d4d7"),
  Common_name = c(Anathera = "#527f9e", Meyeri = "#9c86bf", Strangulata = "#00F5FF",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F")
)
pheatmap(data,cluster_rows = F,cutree_cols = 7,color = cols,fontsize=4.5,annotation_col = annotation_c,
         annotation_colors = ann_colors)


####这是冰川的数据
library(ggplot2)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder")
chr1B = read.table("find_chr1B.txt",header = T)
ggplot(chr1B, aes(x=chr1B$location, y=chr1B$LR)) +
  geom_line()

data <- read.table("LR04stack.txt",header=T,stringsAsFactors = F,sep="\t")
ggplot(data, aes(x=data$Time..ka., y=data$Benthic.d18O..per.mil.)) +
  geom_line() +
  scale_x_log10()+scale_y_reverse() +
  geom_hline(aes(yintercept=4), colour="#990000", linetype="dashed")


#####现在是把Q基因的情况画haplotype sharing
library(RColorBrewer)
display.brewer.all()
cols= brewer.pal(8,"Oranges")
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder")
chr26Q_heatmap = read.table("chr26Q.cp2.chunkcounts.out2",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:1)])
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5)

library(RColorBrewer)
cols= brewer.pal(6,"Oranges")
CpainterD = read.table("unlinkedblockD.unnamed.chunklengths.out",header=T)
data11 = na.omit(CpainterD)
data1 = abs(log(data11[,-(1:1)]))
data = as.matrix(data1)
df <- data[is.finite(rowSums(data)),]
#pheatmap(df,cluster_rows = F,cutree_cols = 3,cluster_rows = 3,kmeans_k = 10,color = cols,fontsize=4.5)
pheatmap(df,cluster_rows = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=4.5)


####
library(heatmap.plus)
library(pheatmap)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/Cpainter/top30")
CpainterD = read.table("ChromoCombine_top30_blockD.chunklengths.out",header=T)
data11 = na.omit(CpainterD)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
info1 = read.delim("annotationD.csv",head =T,sep=",")
info_head = colnames(CpainterD[,-(1:1)])
index=match(info_head,info1[,2])
info2 = info1[index,]
annotation_c=info2[,c(3,4)]
rownames(annotation_c) <- colnames(data)
info_head = rownames(CpainterD[,-(1:1)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c2=info2[,c(3,4)]
rownames(annotation_c2) <- rownames(data)
ann_colors = list(
  Class = c(DD = "#615f60", ABD_orphon = "#a19d9f",Bread_wheat = "#d9d4d7"),
  Common_name = c(Anathera = "#527f9e", Meyeri = "#9c86bf", Strangulata = "#00F5FF",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F")
)
cols= brewer.pal(8,"Oranges")
pheatmap(data,cutree_cols = 3,cluster_rows = 3,kmeans_k = 10,color = cols,fontsize=4.5)
#pheatmap(data,cluster_rows = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=4.5)
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, annotation_row = annotation_c2,
         annotation_colors = ann_colors)
###画D的一部分
CpainterD = read.table("D.chunklengths.txt",header=T)
data11 = na.omit(CpainterD)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, 
         annotation_colors = ann_colors)

######BB
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/Cpainter/top30")
CpainterB = read.table("ChromoCombine_top30_blockD.chunklengths.out",header=T)
data11 = na.omit(CpainterB)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
info1 = read.delim("annotationB.csv",head =T,sep=",")
info_head = colnames(CpainterB[,-(1:1)])
index=match(info_head,info1[,2])
info2 = info1[index,]
annotation_c=info2[,c(3,4)]
rownames(annotation_c) <- colnames(data)
info_head = rownames(CpainterB[,-(1:1)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c2=info2[,c(3,4)]
rownames(annotation_c2) <- rownames(data)
brewer.pal(8,"Dark2")
ann_colors = list(
  Class = c(SS = "#615f60", AABB = "#a19d9f",AABBDD = "#d9d4d7"),
  Common_name = c(Speltoides = "#527f9e", Wild_emmer = "#9c86bf", Domesticated_emmer = "#00F5FF",
                  Ispahanicum = "#1B9E77", Georgian_wheat = "#1B9E77", Rivet_wheat = "#D95F02",
                  Polish_wheat = "#7570B3", Persian_wheat = "#E7298A", Durum = "#66A61E",Khorasan_wheat = "#E6AB02",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F")
)
cols= brewer.pal(8,"Oranges")
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, annotation_row = annotation_c2,
         annotation_colors = ann_colors)
###画D的一部分
CpainterD = read.table("B.chunklengths.txt",header=T)
data11 = na.omit(CpainterD)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, 
         annotation_colors = ann_colors)

######AA
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/Cpainter/top30")
CpainterA = read.table("ChromoCombine_top30_blockA.chunklengths.out",header=T)
data11 = na.omit(CpainterA)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
info1 = read.delim("annotationA.csv",head =T,sep=",")
info_head = colnames(CpainterA[,-(1:1)])
index=match(info_head,info1[,2])
info2 = info1[index,]
annotation_c=info2[,c(3,4)]
rownames(annotation_c) <- colnames(data)
info_head = rownames(CpainterA[,-(1:1)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c2=info2[,c(3,4)]
rownames(annotation_c2) <- rownames(data)
brewer.pal(8,"Dark2")
ann_colors = list(
  Class = c(AA = "#615f60", AABB = "#a19d9f",AABBDD = "#d9d4d7"),
  Common_name = c(Wild_einkorn = "#527f9e", Domesticated_einkorn = "#75968d",Urartu = "#967594",
                  Wild_emmer = "#9c86bf", Domesticated_emmer = "#00F5FF",
                                                                                                                             Ispahanicum = "#1B9E77", Georgian_wheat = "#1B9E77", Rivet_wheat = "#D95F02",
                  Polish_wheat = "#7570B3", Persian_wheat = "#E7298A", Durum = "#66A61E",Khorasan_wheat = "#E6AB02",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F")
)
cols= brewer.pal(8,"Oranges")
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, annotation_row = annotation_c2,
         annotation_colors = ann_colors)
###画D的一部分
CpainterA = read.table("A.chunklengths.txt",header=T)
data11 = na.omit(CpainterD)
data1 = log(data11[,-(1:1)] + 1)
data = as.matrix(data1)
pheatmap(data,cluster_rows = T,cluster_cols = T,cutree_rows = 7,cutree_cols =8,color = cols,fontsize=6,
         annotation_col = annotation_c, 
         annotation_colors = ann_colors)



######现在是两个文件的GO比较
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/GO/B")
all=read.table("Go_B_sort2.txt",sep = "\t",header=F)
A1 = read.table("wild1.txt",sep="\t",header =F)
RR = all[match(all$V1,A1$V1),1]
RR = as.character(RR)
for (i in c(1:length(RR))){
  if (is.na(RR[i])){
    RR[i] <- "0"
  }else
    RR[i] <- "1"
}
write.table(RR,"test.txt",sep="\t",col.names = F,row.names = F,quote=F)

##画图
library(heatmap.plus)
library(pheatmap)
library(RColorBrewer)
##A
display.brewer.all()
brewer.pal(9,'BuPu')
cols = c("#e9d5f0","#88419D")
data <- read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/GO/Gene_Go_A.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data = as.matrix(data)
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=7)
#data11 = t(data)
#pheatmap(data11,cluster_rows = F,cluster_cols = F,color = cols,fontsize=6)
##B
display.brewer.all()
brewer.pal(9,'Greens')
cols = c("#a3d9b8","#006D2C")
data <- read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/GO/Gene_Go_B.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data = as.matrix(data)
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=8)
data11 = t(data)
#pheatmap(data11,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=6)
##D
display.brewer.all()
brewer.pal(9,'YlOrRd')
cols = c("#f7e6d0","#FEB24C")
data <- read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/volcanofinder/GO/Gene_Go_D.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data = as.matrix(data)
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=8)
data11 = t(data)



