source("genomics_general/jackknife.R") 
freq_table = read.table("data/A.tsv.gz", header=T, as.is=T) 
info = read.table("AF_090919.info.txt",head=T) 
info[,1:4] = sapply(info[,1:4],as.character) 
for (i in 1:nrow(info)){ 
  P1 = info[i,2] 
  P2 = info[i,3] 
  P3 = info[i,4] 
  co = paste("A0910",info[i,1],sep="") 
  re = getDstat(freq_table,P1,P2,P3,co) 
  info[i,5] = re$D 
  info[i,6] = re$Z 
  info[i,7] = re$f 
  info[i,8] = re$f_CI_lower 
  info[i,9] = re$f_CI_upper 
} 
write.table(info,"AF_090919.info.all.txt",col.names = T,row.names = F,quote=F,sep="\t") 


###现在想要看的是在wild emmer形成之后，speltoides对dom emmer的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
cols1 = colorRampPalette(c("white","#FF0000"))(30)
cols2 = colorRampPalette(c("white","#fdd662"))(30)
cols4 = colorRampPalette(c("white","#556B2F"))(30)
# colEA1 = cols2[25]   "#FDDD7D"
# colWA1 = cols1[8]   "#FFC1C1"
# colEU1 = cols4[15] "#ACB79A"
colEA = rgb(253/255,221/255,125/255,alpha=0.7) 
colWA = rgb(255/255,193/255,193/255,alpha=0.6) 
colEU = rgb(172/255,183/255,154/255,alpha=1) 
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Speltoides/")
xpxlrA = read.table("fd_Speltoides_domemmer.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,1), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = colEU,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,0.5))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Speltoides to dom emmer (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 2)
  }
}
axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = 400)

##刚刚画的图，305843558-308909914的这个区域的introgression的位置是introgression比较高的位置--3M
##300562860-327838165--20M
##现在是把这些位置画haplotype
library(heatmap.plus)
library(pheatmap)
######这是测试的第一个
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Speltoides")
chr26Q_heatmap = read.table("fd_1B_1M.txt",header=T)
data = as.matrix(chr26Q_heatmap[,-(1:4)])
info1 = read.delim("allB_subspecies.csv",head =T,sep=",")
info_head = colnames(chr26Q_heatmap[,-(1:4)])
index=match(info_head,info1[,1])
info2 = info1[index,]
annotation_c=info2[,c(2,3)]
rownames(annotation_c) <- colnames(data)
ann_colors = list(
  Class = c(SS = "#615f60", AABB = "#a19d9f",AABBDD = "#d9d4d7"),
  Common_name = c(Speltoides = "#527f9e", Rivet_wheat = "#9c86bf", Polish_wheat = "#00F5FF",
                  Wild_emmer = "#1B9E77", Khorasan_wheat = "#1B9E77", Domesticated_emmer = "#D95F02",
                  Persian_wheat = "#7570B3", Ispahanicum = "#E7298A", Georgian_wheat = "#66A61E",Durum = "#E6AB02",
                  Club_wheat ="#FDD662",Macha = "#8F9E76", Spelt ="#FF8C8C",
                  Indian_dwarf_wheat = "#1F78B4", Yunan_wheat = "#BEBADA", Tibetan_semi_wild = "#FB8072",
                  Xinjiang_wheat = "#FFFFB3", Vavilovii = "#D9D9D9", Landrace = "#BC80BD",Cultivar = "#FFED6F",Synthetic = "#d9d4d7")
)
cols = c("#FFF5EB","#FCCDE5","#fa75b9")
pheatmap(data,cluster_rows = F,cutree_cols = 3,color = cols,fontsize=4.5,annotation_col = annotation_c,
         annotation_colors = ann_colors)


library(vioplot)
library(viridis)
cols = viridis(4)
cols1 = colorRampPalette(c("white",cols[1]))(30)
cols2 = colorRampPalette(c("white",cols[2]))(30)
cols4 = colorRampPalette(c("white",cols[4]))(30)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/IntroV1/")
par(mfrow=c(1,2), oma=c(5,2, 2, 2), mar=c(5,4,2,2), cex=1,bty="n")
#Uratu_tatraploids
all = read.table("general_Uratu_tatraploids.txt",head =T)
all1 = all
all$ID = factor(all1$ID,levels = c("fd_Uratu_North_wildemmer","fd_Uratu_Domesticated_emmer","fd_Uratu_Ispahanicum",
                                   "fd_Uratu_Georgian_wheat","fd_Uratu_Rivet_wheat","fd_Uratu_Polish_wheat",
                                   "fd_Uratu_Persian_wheat","fd_Uratu_Durum","fd_Uratu_Khorasan_wheat"))
aa = boxplot(all$fd~all$ID,col=c(rep(cols1[20],4),rep(cols2[30],5)),cex.axis=0.0001,yaxt ="n",outline=F,ylim=c(0,0.25))
axis(2)
mtext("fd",side = 2 ,las = 0, line  = 3 ,cex = 1.3,at = 0.1)
mtext("Uratu_tatraploids",side = 3 ,las = 0, line  = 0 ,cex = 1.1,at = 0.6)
#Speltoides_tatraploids
all = read.table("general_Speltoides_tatraploids.txt",head =T)
all1 = all
all$ID = factor(all1$ID,levels = c("fd_Speltoides_North_wildemmer","fd_Speltoides_Domesticated_emmer","fd_Speltoides_Ispahanicum",
                                   "fd_Speltoides_Georgian_wheat","fd_Speltoides_Rivet_wheat","fd_Speltoides_Polish_wheat",
                                   "fd_Speltoides_Persian_wheat","fd_Speltoides_Durum","fd_Speltoides_Khorasan_wheat"))
aa = boxplot(all$fd~all$ID,col=c(rep(cols1[20],4),rep(cols2[30],5)),cex.axis=0.0001,yaxt ="n",outline=F,ylim=c(0,0.25))
axis(2)
mtext("fd",side = 2 ,las = 0, line  = 3 ,cex = 1.3,at = 0.1)
mtext("Speltoides_tatraploids",side = 3 ,las = 0, line  = 0 ,cex = 1.1,at = 0.6)

###proportion
pro1 = read.table("proportion_Uratu_tatraploids.txt",head =T)
plot(pro1$ratio,type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,ylim = c(0,0.04),axes = F)
axis(4,cex.axis = 1,col="#D4005F")
pro2 = read.table("proportion_Speltoides_tatraploids.txt",head =T)
plot(pro2$ratio+0.01,type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,ylim = c(0,0.04),axes = F)
axis(4,cex.axis = 1,col="#D4005F")

#######################################
par(mfrow=c(1,2), oma=c(5,2, 2, 2), mar=c(5,4,2,2), cex=1,bty="n")
#Tatraploids_Hexa
all = read.table("general_Tatraploids_Hexa.txt",head =T)
all1 = all
all$ID = factor(all1$ID,levels = c("fd_Tatraploids_Yunan_wheat","fd_Tatraploids_Tibetan_semi_wild","fd_Tatraploids_Xinjiang_wheat",
                                   "fd_Tatraploids_Club_wheat","fd_Tatraploids_Spelt","fd_Tatraploids_Macha",
                                   "fd_Tatraploids_Vavilovii","fd_Tatraploids_Landrace","fd_Tatraploids_Cultivar"))
aa = boxplot(all$fd~all$ID,col=c(rep(cols1[7],3),rep(cols2[17],4),rep(cols4[10],2)),cex.axis=0.0001,yaxt ="n",outline=F,ylim=c(0,1))
axis(2)
mtext("fd",side = 2 ,las = 0, line  = 3 ,cex = 1.3,at = 0.5)
mtext("Tatraploids_Hexa",side = 3 ,las = 0, line  = 0 ,cex = 1.1,at = 0.6)
#Strangulata_Hexa
all = read.table("general_Strangulata_Hexa.txt",head =T)
all1 = all
all$ID = factor(all1$ID,levels = c("fd_Strangulata_Yunan_wheat","fd_Strangulata_Tibetan_semi_wild","fd_Strangulata_Xinjiang_wheat",
                                   "fd_Strangulata_Club_wheat","fd_Strangulata_Spelt","fd_Strangulata_Macha",
                                   "fd_Strangulata_Vavilovii","fd_Strangulata_Landrace","fd_Strangulata_Cultivar"))
aa = boxplot(all$fd~all$ID,col=c(rep(cols1[7],3),rep(cols2[17],4),rep(cols4[10],2)),cex.axis=0.0001,yaxt ="n",outline=F,ylim=c(0,1))
axis(2)
mtext("fd",side = 2 ,las = 0, line  = 3 ,cex = 1.3,at = 0.1)
mtext("Strangulata_Hexa",side = 3 ,las = 0, line  = 0 ,cex = 1.1,at = 0.6)

###proportion
pro1 = read.table("proportion_Tatraploids_Hexa.txt",head =T)
plot(pro1$ratio,type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,ylim = c(0,0.25),axes = F)
axis(4,cex.axis = 1,col="#D4005F")
pro2 = read.table("proportion_Strangulata_Hexa.txt",head =T)
plot(pro2$ratio+0.01,type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,ylim = c(0,0.25),axes = F)
axis(4,cex.axis = 1,col="#D4005F")







