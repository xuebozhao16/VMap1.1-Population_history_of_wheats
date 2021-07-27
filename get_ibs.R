setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/IBSdis_fst/testcode/")
group = read.table("wildemmer_domemmer_20.txt",header = F)
group1 = group[1:10,]
group2 = group[11:20,]
info = read.table("chr1_1M.bed",head=F)
all_ibs = data.frame(chr=integer(),pos1 = integer(),pos2 = integer(), ibs_value= numeric())
for(i in 1:length(info$V1)){
  print(paste("Reading file ",info[i,1:3],sep=""))
  name = paste(info[i,1],".",info[i,2],"-",info[i,3],".ibs.txt",sep="")
  ibs = read.table(name,header = T)
  index1 = match(group1,ibs$Dxy)
  ibs_rm_group2 = ibs[,c(1,index1+1)]
  index2 = match(group2,ibs_rm_group2$Dxy)
  ibs_rm_group1 = ibs_rm_group2[index2,]
  ibs_frame = as.data.frame(ibs_rm_group1[,-1])
  meanibs = sum(ibs_frame)/100
  all_ibs[i,1] = info[i,1]
  all_ibs[i,2] = info[i,2]
  all_ibs[i,3] = info[i,3]
  all_ibs[i,4] = meanibs
}
write.table(all_ibs,"./wildemmer_domemmer_ibs.txt",col.names = T,row.names = F,quote=F,sep="\t")

#####开始整理IBS的图
#####wildemmer-domemmer的群体的IBS
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
colEA = rgb(253/255,221/255,125/255,alpha=0.7) 
colWA = rgb(255/255,193/255,193/255,alpha=0.9) 
colEU = rgb(172/255,183/255,154/255,alpha=1) 
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/IBSdis_fst/pop_ibs")
xpxlrA = read.table("Chr_wildemmer_domemmer_ibs.txt",head =F)
xpxlrB = read.table("Chr_domemmer_free_ibs.txt",head =F)
xpxlrC = read.table("Chr_Free_Land_ibs.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = "blue",type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,0.3))
  score2 = xpxlrB[xpxlrB[,1] == c,3]
  pos2 = xpxlrB[xpxlrB[,1] == c,2]
  lines(pos2/1000000,score2,col = "orange")
  score3 = xpxlrC[xpxlrC[,1] == c,3]
  pos3 = xpxlrC[xpxlrC[,1] == c,2]
  lines(pos3/1000000,score3,col = "green")
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col=col)
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("IBS distence",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 1.2)
  }
}
axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.6 ,cex = 1,at = -50)

####整理box-plot
library(RColorBrewer)
par(mfrow=c(3,1))
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/IBSdis_fst/pop_ibs")
wildemmer_domemmer = read.table("Chr_wildemmer_domemmer_ibs.txt",head =F)
compare1 = data.frame(wildemmer_domemmer[,c(1,3)])
boxplot(V3~V1,compare1,col= brewer.pal(7,'Blues'),outline = F,xlab = "CHR",ylab = "IBS distence",
        ylim = c(0,0.35), main="Wildemmer vs Domemmer")
domemmer_free = read.table("Chr_domemmer_free_ibs.txt",head =F)
compare2 = data.frame(domemmer_free[,c(1,3)])
boxplot(V3~V1,compare2,col= brewer.pal(7,'Oranges'),outline = F,xlab = "CHR",ylab = "IBS distence",
        ylim = c(0,0.35), main="Domemmer vs Free-threshing")
Free_Land = read.table("Chr_Free_Land_ibs.txt",head =F)
compare3 = data.frame(Free_Land[,c(1,3)])
boxplot(V3~V1,compare3,col= brewer.pal(7,'Greens'),outline = F,xlab = "CHR",ylab = "IBS distence",
        ylim = c(0,0.35), main="Free-threshing vs Landrace AB")

par(mfrow=c(2,1))
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/IBSdis_fst/pop_ibs")
wildemmer_domemmer = read.table("Chr_wildemmer_domemmer_ibs.txt",head =F)
compare1 = data.frame(wildemmer_domemmer[,c(1,3)])
boxplot(V3~V1,compare1,col= brewer.pal(7,'Blues'),outline = T,xlab = "CHR",ylab = "IBS distence",
        ylim = c(0,0.7), main="Wildemmer vs Domemmer")
domemmer_free = read.table("Chr_domemmer_free_ibs.txt",head =F)
compare2 = data.frame(domemmer_free[,c(1,3)])
boxplot(V3~V1,compare2,col= brewer.pal(7,'Oranges'),outline = T,xlab = "CHR",ylab = "IBS distence",
        ylim = c(0,0.7), main="Domemmer vs Free-threshing")

#####25&-75%
par(mfrow=c(1,3))
plot(density(wildemmer_domemmer$V3), xlim = c(0,0.6),xlab = "IBS distence", main="Wildemmer vs Domemmer")
summary(wildemmer_domemmer$V3)
abline(v = 0.0816, col = "gray60")
abline(v = 0.1289, col = "gray60")
plot(density(domemmer_free$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Domemmer vs Free-threshing")
summary(domemmer_free$V3)
abline(v = 0.05226, col = "gray60")
abline(v = 0.11703, col = "gray60")
plot(density(Free_Land$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Free-threshing vs Landrace AB")
summary(Free_Land$V3)
abline(v = 0.03784, col = "gray60")
abline(v = 0.11081, col = "gray60")
#####35&-65%
par(mfrow=c(1,3))
plot(density(wildemmer_domemmer$V3), xlim = c(0,0.6),xlab = "IBS distence", main="Wildemmer vs Domemmer")
sort_data = sort(wildemmer_domemmer$V3)
aa = sort_data[ceiling(length(wildemmer_domemmer$V3)*0.35)]
bb = sort_data[ceiling(length(wildemmer_domemmer$V3)*0.65)]
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")
plot(density(domemmer_free$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Domemmer vs Free-threshing")
sort_data = sort(domemmer_free$V3)
aa = sort_data[ceiling(length(domemmer_free$V3)*0.35)]
bb = sort_data[ceiling(length(domemmer_free$V3)*0.65)]
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")
plot(density(Free_Land$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Free-threshing vs Landrace AB")
sort_data = sort(Free_Land$V3)
aa = sort_data[ceiling(length(Free_Land$V3)*0.35)]
bb = sort_data[ceiling(length(Free_Land$V3)*0.65)]
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")
#####peak
par(mfrow=c(1,3))
plot(density(wildemmer_domemmer$V3), xlim = c(0,0.6),xlab = "IBS distence", main="Wildemmer vs Domemmer")
max = which.max(density(wildemmer_domemmer$V3)$y)
peak=density(wildemmer_domemmer$V3)$x[max]
peak
sort_data = sort(wildemmer_domemmer$V3)
max_locate = which(sort_data == 0.1078)
aa = sort_data[max_locate-ceiling(length(wildemmer_domemmer$V3)*0.15)]
bb = sort_data[max_locate+ceiling(length(wildemmer_domemmer$V3)*0.15)]
abline(v = peak, col = "red")
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")
plot(density(domemmer_free$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Domemmer vs Free-threshing")
max = which.max(density(domemmer_free$V3)$y)
peak=density(domemmer_free$V3)$x[max]
peak
sort_data = sort(domemmer_free$V3)
max_locate = which(sort_data == 0.0926)
aa = sort_data[max_locate-ceiling(length(domemmer_free$V3)*0.15)]
bb = sort_data[max_locate+ceiling(length(domemmer_free$V3)*0.15)]
abline(v = peak, col = "red")
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")
plot(density(Free_Land$V3), xlim = c(0,0.6),xlab = "IBS distence",main="Free-threshing vs Landrace AB")
max = which.max(density(Free_Land$V3)$y)
peak=density(Free_Land$V3)$x[max]
peak
sort_data = sort(Free_Land$V3)
max_locate = which(sort_data == 0.0845)
aa = sort_data[max_locate-ceiling(length(Free_Land$V3)*0.15)]
bb = sort_data[max_locate+ceiling(length(Free_Land$V3)*0.15)]
abline(v = peak, col = "red")
abline(v = aa, col = "gray60")
abline(v = bb, col = "gray60")



