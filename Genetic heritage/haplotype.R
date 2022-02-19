setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/permutation/Crispr_gene/TraesCS5A02G233700/SNP")
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
#install.packages("ggseqlogo")
library(ggplot2)
library(ggseqlogo)
data(ggseqlogo_sample)
fasta = "ACGTATG
ATGTATG
ACGTATG
ACATATG
ACGTACG"
ggplot()+geom_logo(seqs_dna$MA0001.1)+theme_logo()

##OTUB1表达量的情况
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1")
par(mfrow=c(3,1),oma=c(4.5,4, 4, 2.5), mar=c(2,1,1,1), cex=1)
layout(matrix(c(1:3), 3, 1, byrow = TRUE), widths=rep(1,3), heights=c(rep(1,3),0.002))
gene7A=read.delim("expression_7A.csv",head=T,sep = ",")
data1 = gene7A[,2]
data2= gene7A[,3]
aa=barplot(as.matrix(t(data1)),ylab = "7A expression Level(TPM)",col = "lavender",ylim = c(0,30))
arrows(aa,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7A expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 15)
gene7B=read.delim("expression_7B.csv",head=T,sep = ",")
data1 = gene7B[,2]
data2= gene7B[,3]
bb=barplot(as.matrix(t(data1)),ylab = "7B expression Level(TPM)",col = "DarkSeaGreen1",ylim = c(0,30))
arrows(bb,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7B expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 15)
gene7D=read.delim("expression_7D.csv",head=T,sep = ",")
data1 = gene7D[,2]
data2= gene7D[,3]
dd=barplot(as.matrix(t(data1)),ylab = "7D expression Level(TPM)",col = "LightGoldenrod3",ylim = c(0,30))
arrows(dd,data1-data2,aa,data1+data2,length=0.05, angle=90, code=3)
mtext("7D expression Level(TPM)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 15)
cap1 = gene7D[,1]
text(dd, y=par()$usr[3]-2,srt=30, xpd=T, adj=1, labels = cap1,cex = 1)

###OTUB1单倍型图谱
cols = c("#FFF5EB","#80B1D3","#FCCDE5")
library(heatmap.plus)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/cooperation_fu/OTUB1")
chr26Q_heatmap = read.table("chr37_TraesCS7A02G263900_2K.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
chr26Q_heatmap = read.table("chr39_TraesCS7B02G161900_2K.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)
chr26Q_heatmap = read.table("chr41_TraesCS7D02G264800_2K.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
heatmap.plus(t(data),na.rm = T,scale = "none",col=cols,Colv = NA,Rowv = NA, cexRow=0.4,cexCol = 0.4)












