cols1 = "#8da0cb"
cols2 = "#00CDCD"

cols3 = "#e78ac3"
cols4 = "#CD919E"

cols5 = "#66c2a5"
cols6 = "#6E8B3D"

cols7 = "#fc8d62"
cols8 = "#FFA500"
###Dom emmer对spelt的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Spelt")
xpxlrA = read.table("fd_domemmer_spelt.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols1,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Domesticated emmer to spelt (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)

###Landrace对spelt的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Spelt")
xpxlrA = read.table("fd_landrace_spelt.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols2,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Landrace to spelt (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)


###Dom emmer对Macha的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Macha")
xpxlrA = read.table("Smooth2M_domemmer_Macha.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols3,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Domesticated emmer to Maccha (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)

###Landrace对Macha的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Macha")
xpxlrA = read.table("fd_landrace_Macha.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols4,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Landrace to Macha (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)


###Polish_wheat对emmer对Xinjiang_wheat的渗入的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Xinjiang")
xpxlrA = read.table("Smooth2M_Polish_wheat_Xinjiang_wheat.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols5,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Polish_wheat to Xinjiang_wheat (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)

###Landrace对emmer对Xinjiang_wheat的渗入的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Xinjiang")
xpxlrA = read.table("Smooth2M_landrace_Xinjiang_wheat.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols6,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("landrace to Xinjiang_wheat (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)


###Rivet_wheat对Persian的渗入的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Persian")
xpxlrA = read.table("Smooth2M_Rivet_wheat_Persian_wheat.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols7,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Smooth2M_Rivet_wheat to Persian_wheat.txt (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)


###Landracec对Persian的渗入的渗入
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/introgression/Persian")
xpxlrA = read.table("Smooth2M_landrace_Persian_wheat.txt",head =F)
max(xpxlrA$V3)
chr = xpxlrA[,1]
chr.num1 <- unique(chr)
chr.num=sort(chr.num1)
par(mfrow=c(7,2), oma=c(4.5,3, 0, 0), mar=c(0.3,1.6,0.3,1), cex=1)
for(c in chr.num){
  score = xpxlrA[xpxlrA[,1] == c,3]
  pos = xpxlrA[xpxlrA[,1] == c,2]
  plot(pos/1000000,score,axes = F,col = cols8,type="l",lwd=1.5,xlim=c(0,800),ylim = c(0,1))
  axis(2,cex.axis =0.8)
  index = cent$chr==c
  points(cent[index,2],0.0002,pch=20,cex=2,col="grey")
  mtext(c,side = 3 ,las = 1, line  = -1.5 ,cex = 1,at = 380)
  if(c == "chr7A"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    #mtext("(Mb)",side = 1 ,las = 1, line  = 2.1 ,cex = 0.8,at = 830)
    mtext("Landrace to Persian_wheat.txt (fd)",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 4)
  }
  if(c == "chr7B"){
    axis(1, pos = -0.1,cex.axis =0.8,at = c(0, 200, 400, 600,800))
  }
}
mtext("Genomic position (Mb)",side = 1 ,las = 1, line  = 2.5 ,cex = 1,at = -100)



















