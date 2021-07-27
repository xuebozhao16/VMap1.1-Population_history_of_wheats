##gene tree的代码整理
##统计每个基因里面有多少个SNP，画一个density plot
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree")
SNPcount = read.table("chrall_genesnpcount.txt",head =F)
hist(SNPcount$V2, breaks=1000, col="red",xlim = c(0,100),ylab="countNum",xlab="Snp in per gene",main = "")
mean(SNPcount$V2)
sd(SNPcount$V2)
d = density(SNPcount$V2)
plot(d)
plot(d,xlim = c(0,50))

##
#library(pegas)
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
h <- haplotype(x)
net <- haploNet(h)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8)

setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/genotype_network")
x = read.dna(file = "Bexontree_withBarleyNoS.fasta" ,format = "fasta")
y <- x[sample(392, size =10, replace = TRUE), ]
h = haplotype(y)
net = haploNet(h)
#plot(net)
plot(net, size=attr(net, "freq"), scale.ratio = 100, cex = 0.1,labels  = T)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 5)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 1)
plot(net, size=attr(net, "freq"), scale.ratio = 6, cex = 1)
plot(net, size=attr(net, "freq"), scale.ratio = 6, cex = 2)
table(rownames(y))
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=rownames(y)[values])
)
ind.hap[1:1, 1:60]
plot(net, size=attr(net, "freq"), scale.ratio = 40 , cex = 0.1, pie=ind.hap,show.mutation =1)
plot(net, size=attr(net, "freq"), scale.ratio = 40 , cex = 0.1, pie=ind.hap,show.mutation =6,labels = TRUE,col.link = "black",threshold = c(1, 10))
legend(30,2, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20,cex = 0.5)
wrong.pop<-rep(letters[1:272], each=22)
ind.hap2<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=wrong.pop[values])
)
plot(net, size=attr(net, "freq"), scale.ratio = 40, cex = 0.1, pie=ind.hap2)
legend(50,20, colnames(ind.hap2), col=rainbow(ncol(ind.hap2)), pch=20)

##depth 和 sd 的分布图
# mean(SS$V1)
# dep
# dep[,1:ncol(dep)] = sapply(dep[,1:ncol(dep)],as.character)
# dep[,1:ncol(dep)] = sapply(dep[,1:ncol(dep)],as.numeric)
# apply(dep[,-(1:2)], 2, sd)
# apply(dep[,-(1:2)], 2, mean)

setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/bwamap/depth_SD")
##现在画的是row bam文件的depth
## pdf(file="DepthandSD.pdf",height=10,width=10)
SS = read.table("all_100k.depth",header = F)
red =  rgb(1,0,0,0.05)
plot(SS,col =red,pch=19,xlab = "Depth",ylab = "sd",xlim = c(0,15),ylim = c(0,10),main="S genome depth")
d = density(SS$V1)
plot(d,col =red,xlim = c(0,30))
hist(SS$V1,col =red, breaks=10000,xlim = c(0,30))
##现在是从row的VCF里面找到depth
SS = read.table("depthSD_chr_all_VCF.txt",header = F)
red =  rgb(1,0,0,0.005)
plot(SS,col =red,pch=19,xlab = "Depth",ylab = "sd",xlim = c(0,50),ylim = c(0,20))
plot(density(SS$V1),col ="red",xlim = c(0,20),main="Depth density",xlab = "Depth")
hist(SS$V1,col =red, breaks=10000,xlim = c(0,30),main="Main depth",xlab = "Depth")
plot(density(SS$V2),col ="red",xlim = c(0,30),main="SD density",xlab = "SD")
hist(SS$V2,col =red, breaks=10000,xlim = c(0,20),main="SD depth",xlab = "SD")
#现在用500K的再试一下
SS = read.table("depthSD_chr_all_VCF_500K.txt",header = F)
red =  rgb(1,0,0,0.01)
plot(SS,col =red,pch=19,xlab = "Depth",ylab = "sd",xlim = c(0,50),ylim = c(0,30))
plot(density(SS$V1),col ="red",xlim = c(0,20),main="Depth density",xlab = "Depth")
hist(SS$V1,col =red, breaks=10000,xlim = c(0,30),main="S genome raw SNP main depth",xlab = "Depth")
plot(density(SS$V2),col ="red",xlim = c(0,30),main="SD density",xlab = "SD")
hist(SS$V2,col =red, breaks=10000,xlim = c(0,20),main="SD depth",xlab = "SD")

###现在是画E3的maf的分布图
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/checkD/E3_maf")
par(mfrow=c(1,3), oma=c(5,6, 1, 1), mar=c(5,5,1,1), cex=1)
E3_A = read.table("E3_100K_Alineage_siteQCfile.txt",header = F)
hist(E3_A$V5,breaks = 100,col = "red",xlab = "MAF",main = "A lineage")
E3_B = read.table("E3_100K_Blineage_siteQCfile.txt",header = F)
hist(E3_B$V5,breaks = 100,col = "red",xlab = "MAF",main = "B lineage")
E3_D = read.table("E3_100K_Dlineage_siteQCfile.txt",header = F)
hist(E3_D$V5,breaks = 100,col = "red",xlab = "MAF",main = "D lineage")

######A
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/checkD/E3_maf/A")
par(mfrow=c(2,4), oma=c(5,6, 1, 1), mar=c(5,5,1,1), cex=1)
hist(E3_A$V5,breaks = 100,col = "red",xlab = "MAF",main = "A lineage")
E3_A1 = read.table("A111_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A111")
E3_A1 = read.table("A110_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A110")
E3_A1 = read.table("A101_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A101")
E3_A1 = read.table("A100_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A100")
E3_A1 = read.table("A011_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A011")
E3_A1 = read.table("A010_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A010")
E3_A1 = read.table("A001_site10K.txt",header = F)
hist(E3_A1$V5,breaks = 100,col = "red",xlab = "MAF",main = "A001")
######B
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/checkD/E3_maf/B")
par(mfrow=c(2,4), oma=c(5,6, 1, 1), mar=c(5,5,1,1), cex=1)
hist(E3_B$V5,breaks = 100,col = "red",xlab = "MAF",main = "B lineage")
E3_B1 = read.table("B111_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B111")
E3_B1 = read.table("B110_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B110")
E3_B1 = read.table("B101_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B101")
E3_B1 = read.table("B100_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B100")
E3_B1 = read.table("B011_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B011")
E3_B1 = read.table("B010_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B010")
E3_B1 = read.table("B001_site10K.txt",header = F)
hist(E3_B1$V5,breaks = 100,col = "red",xlab = "MAF",main = "B001")
######D
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/checkD/E3_maf/D")
par(mfrow=c(2,4), oma=c(5,6, 1, 1), mar=c(5,5,1,1), cex=1)
hist(E3_D$V5,breaks = 100,col = "red",xlab = "MAF",main = "D lineage")
E3_D1 = read.table("D111_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D111")
E3_D1 = read.table("D110_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D110")
E3_D1 = read.table("D101_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D101")
E3_D1 = read.table("D100_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D100")
E3_D1 = read.table("D011_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D011")
E3_D1 = read.table("D010_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D010")
E3_D1 = read.table("D001_site10K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "D001")

######checkD and V11
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/checkD/E3_maf/D")
par(mfrow=c(1,2), oma=c(5,6, 1, 1), mar=c(5,5,1,1), cex=1)
hist(E3_D$V5,breaks = 100,col = "red",xlab = "MAF",main = "Vmap1.1 D lineage")
E3_D1 = read.table("V11_checkD_siteQCfile_100K.txt",header = F)
hist(E3_D1$V5,breaks = 100,col = "red",xlab = "MAF",main = "Vmap1.0 D lineage")

















