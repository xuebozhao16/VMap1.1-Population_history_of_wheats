source('~/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/DensityPlot_snpdensi.R')
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity")
par(mfrow=c(2,2), oma=c(3,3, 0.5, 1), mar=c(2,2,2,2), cex=1)
cent = read.table("/Users/xuebozhao/Document/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
S = read.table("A.den",head =T)
head(S)
max(S[,4])
min(S[,4])
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'Blues')
density(S[,3])
DensityPlot_snpdensi(S,cent,legend.min = 10,legend.max = 600,legend.len = 8,main = "A Lineage SNP density",col = col1,legend.x.intersp=0.6)
S = read.table("B.den",head =T)
DensityPlot_snpdensi(S,cent,legend.min = 10,legend.max = 600,legend.len = 8,main = "B Lineage SNP density",col = col1,legend.x.intersp=0.6)
S = read.table("D.den",head =T)
DensityPlot_snpdensi(S,cent,legend.min = 10,legend.max = 600,legend.len = 8,main = "D Lineage SNP density",col = col1,legend.x.intersp=0.6)

####现在是为了得到depth分布的QQ plot
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/depth")
par(mfrow=c(1,2))
count=read.table("snpcount.txt",header = F)
s1 = sort(count$V4)
s2 = sort(count$V5)
plot(s1,s2,xlab = "Vmap1 SNP count",ylab = "Vmap1gentic SNP count")
s3 = sort(count$V6)
s4 = sort(count$V7)
plot(s3,s4,xlab = "fd SNP count",ylab = "fdgentic SNP count")


####现在是为了画4条线
#cat *.snpden > all.snpden
#sort -k1,1n -k2,2n all.snpden > withOutgroup_noSyn.snpden
#sort -k1,1n -k2,2n all.snpden > withOutgroup_Syn.snpden
#sort -k1,1n -k2,2n all.snpden > NoOutgroup_Syn.snpden
#sort -k1,1n -k2,2n all.snpden > NoOutgroup_noSyn.snpden
dev.off()
par(mfrow=c(2,2))
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/4line_compare")
No1 = read.table("NoOutgroup_noSyn.snpden",header = T)
No2 = read.table("NoOutgroup_Syn.snpden",header = T)
with1 = read.table("withOutgroup_noSyn.snpden",header = T)
with2 = read.table("withOutgroup_Syn.snpden",header = T)
hist(No1$SNP_COUNT,breaks = 100,col="red",main = "NoOutgroup_noSyn",xlab = "SNP count")
hist(No2$SNP_COUNT,breaks = 100,col="red",main = "NoOutgroup_Syn",xlab = "SNP count")
hist(with1$SNP_COUNT,breaks = 100,col="red",main = "withOutgroup_noSyn",xlab = "SNP count")
hist(with2$SNP_COUNT,breaks = 100,col="red",main = "withOutgroup_Syn",xlab = "SNP count")
##
source('~/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/DensityPlot_snpdensi.R')
par(mfrow=c(2,2), oma=c(3,3, 0.5, 1), mar=c(2,2,2,2), cex=1)
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
S = read.table("Chr_NoOutgroup_noSyn.snpden",head =T)
head(S)
max(S[,4])
min(S[,4])
library(RColorBrewer)
brewer.pal(9,'Blues')
density(S[,3])
DensityPlot_snpdensi(S,cent,legend.min = 100,legend.max = 30000,legend.len = 10,main = "NoOutgroup_noSyn",col = col1,legend.x.intersp=0.6)
S = read.table("Chr_NoOutgroup_Syn.snpden",head =T)
DensityPlot_snpdensi(S,cent,legend.min = 100,legend.max = 30000,legend.len = 10,main = "NoOutgroup_Syn",col = col1,legend.x.intersp=0.6)
S = read.table("Chr_withOutgroup_noSyn.snpden",head =T)
DensityPlot_snpdensi(S,cent,legend.min = 10,legend.max = 1000,legend.len = 10,main = "withOutgroup_noSyn",col = col1,legend.x.intersp=0.6)
S = read.table("Chr_withOutgroup_Syn.snpden",head =T)
DensityPlot_snpdensi(S,cent,legend.min = 10,legend.max = 1000,legend.len = 10,main = "withOutgroup_Syn",col = col1,legend.x.intersp=0.6)


##展示整个基因组的density情况
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
col = "#696969"
cols1 = colorRampPalette(c("white","#FF0000"))(30)
cols2 = colorRampPalette(c("white","#fdd662"))(30)
cols4 = colorRampPalette(c("white","#556B2F"))(30)
colA = rgb(253/255,221/255,125/255,alpha=1) 
colB = rgb(255/255,193/255,193/255,alpha=1) 
colC = rgb(172/255,183/255,154/255,alpha=1) 
##这是裸粒四倍体向landrace渗入
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/2line")
xpxlrA = read.table("Chr_NoOutgroup_noSyn.snpden",header = T)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/2line")
xpxlrB = read.table("Chr_NoOutgroup_Syn.snpden",head =T)
max(xpxlrA$SNP_COUNT)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,3), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = colA,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,45000))
  #plot(pos/1000000,score,axes = F,col = colEU,pch=20,cex=0.7,xlim=c(0,800),ylim = c(0,1))
  score2 = xpxlrB[xpxlrB[,1] == c,3]
  pos2 = xpxlrB[xpxlrB[,1] == c,2]
  lines(pos2/1000000,score2,col = colC)
  #points(pos2/1000000,score2,col = colEA,pch=20,cex=0.7)
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col=col)
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  # if(c == "chr1A"){
  #   legend(570,0.008,c("wild emmer", "domesticated emmer","durum"),col=c("#FDB462", "#B3DE69","#CDB38B"),cex=0.7,lty=1,bty="n",lwd=3)
  # }
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("SNP count",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 180000)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 1,at = -500)
###
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/snp_calling/snpdensity/2line")
dat1 = read.table("a.snpden",head=F)
dat2 = read.table("b.snpden",head=F)
dat3 = read.table("c.snpden",head=F)
dat4 = read.table("d.snpden",head=F)

ia = paste(dat1$V1,dat1$V2,sep="_")
ib = paste(dat2$V1,dat2$V2,sep="_")
ic = paste(dat3$V1,dat3$V2,sep="_")
id = paste(dat4$V1,dat4$V2,sep="_")

as = intersect(ia,ib)
idx1 = match(as,ia)
idx2 = match(as,ib)
plot(dat1[na.omit(idx1),3],dat2[na.omit(idx2),3],pch=19,cex=0.5,xlab="Before synteny filter",ylab = "After synteny fitler")
cor(dat1[na.omit(idx1),3],dat2[na.omit(idx2),3])

colr = rgb(1,0,0,0.2)
as = intersect(id,ic)
idx1 = match(as,id)
idx2 = match(as,ic)
plot(dat4[na.omit(idx1),3],dat3[na.omit(idx2),3],col="black",pch=19,cex=0.5,xlab="Before synteny filter",ylab = "After synteny fitler")
d1 = dat4[na.omit(idx2),]
id1 = paste(d1[,1],d1[,2],sep = "_")
idx = id1==ic

A = c(1,2,7,8,13,14,19,20,25,26,31,32,37,38)
B = c(3,4,9,10,15,16,21,22,27,28,33,34,39,40)
D = c(5,6,11,12,17,18,23,24,29,30,35,36,41,42)

dat5 = dat3[dat3$V1 %in% A,]
dat6 = dat4[dat4$V1 %in% A,]
ie = paste(dat5$V1,dat5$V2,sep="_")
ig = paste(dat6$V1,dat6$V2,sep="_")
as= intersect(ie,ig)
idx1 = match(as,ig)
idx2 = match(as,ie)
points(dat6[na.omit(idx1),3],dat5[na.omit(idx2),3],pch=19,cex=0.5,col=colr)

datA1 = dat3[(dat3$V1 %in% A),]
datA2 = dat4[(dat4$V1 %in% A),]

datB1 = dat3[(dat3$V1 %in% B),]
datB2 = dat4[(dat4$V1 %in% B),]

datD1 = dat3[(dat3$V1 %in% D),]
datD2 = dat4[(dat4$V1 %in% D),]


iA1 = paste(datA1$V1,datA1$V2,sep="_")
iA2 = paste(datA2$V1,datA2$V2,sep="_")

iB1 = paste(datB1$V1,datB1$V2,sep="_")
iB2 = paste(datB2$V1,datB2$V2,sep="_")

iD1 = paste(datD1$V1,datD1$V2,sep="_")
iD2 = paste(datD2$V1,datD2$V2,sep="_")

asA = intersect(iA1,iA2)
asB = intersect(iB1,iB2)
asD = intersect(iD1,iD2)

idA1 = match(asA,iA1)
idA2 = match(asA,iA2)

idB1 = match(asB,iB1)
idB2 = match(asB,iB2)

idD1 = match(asD,iD1)
idD2 = match(asD,iD2)
#
cols = viridis(4)
colr = rgb(68/255,1/255, 84/255,0.5)
colg = rgb(53/255,183/255, 121/255,0.5)
colb= rgb(204/255,186/255, 31/255,0.5)
dev.off()
par(mfrow=c(1,3))
plot(datA2[na.omit(idA2),3],datA1[na.omit(idA1),3],ylim=c(0,22000),col=colr,pch=19,cex=0.5,xlab="Before synteny filter",ylab = "After synteny fitler")
plot(datB2[na.omit(idB2),3],datB1[na.omit(idB1),3],ylim=c(0,22000),pch=19,cex=0.5,col=colg,xlab="Before synteny filter",ylab = "After synteny fitler")
plot(datD2[na.omit(idD2),3],datD1[na.omit(idD1),3],ylim=c(0,22000),pch=19,cex=0.5,col=colb,xlab="Before synteny filter",ylab = "After synteny fitler")

#abline(lm(datD1[na.omit(idD1),3] ~ datD2[na.omit(idD2),3]),col="black")
#lenround(cor(datD2[na.omit(idD2),3],datD1[na.omit(idD1),3]),2)
















