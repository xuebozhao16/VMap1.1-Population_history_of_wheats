library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR")
data1 <- read.table("South_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#data2 <- data1[-which(data1$p == "NaN"),]
data1$SNP <- c(1:330025)
#gwasR <- data2[,c(2,3,4,7)]
gwasR <- data1[,c(7,1,2,5)]
#gwasR$P <- 10^(-gwasR$MeanY)
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
highsnp <- c(213040,213071,213159,213222,213259,213516,213558,213599,213647,213664,213731,213737,213748,213837,213841,213876,213882,213902,213951,213970,213979,213993,214021,214034,214054,214061,214073,214079,214115,214149,214166,214173,214206,214211,329310,329469,329542,329652,329989)
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
# 2） 计算每条chr的初始位置
#chr_pos <- chr_len  %>% 
#  mutate(total = cumsum(chr_len) - chr_len) %>%
#  select(-chr_len) 
#3)计算累计SNP的位置
#Snp_pos <- chr_pos %>%
#  left_join(gwasR, ., by="CHR") %>%
#  arrange(CHR, BP) %>%
#  mutate( BPcum = BP + total)

#查看转化后的数据
#head(Snp_pos,2)

#数据准备完成，开始绘图。
#ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
#  geom_point( aes(color=as.factor(CHR)))
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#p <- ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
#设置点的大小，透明度
#  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
#设置颜色
#  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
#设定X轴
#  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
#去除绘图区和X轴之间的gap
#  scale_y_continuous(expand = c(0, 0) ) +  
#添加阈值线
#  geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size = 1.2, linetype = c("dotted", "twodash")) + 
#设置主题
#  theme_bw() +
#  theme(
#    legend.position="none",
#    panel.border = element_blank(),
#    axis.line.y = element_line(),
#    panel.grid.major.x = element_blank(),
#    panel.grid.minor.x = element_blank()
#  )
#head(gwasR)
#as.data.frame(table(gwasR$CHR))
#manhattan(gwasR)
#gwasR[gwasR$CHR=="1",2] <- "1A"
#gwasR[gwasR$CHR=="2",2] <- "1B"
#gwasR[gwasR$CHR=="3",2] <- "1D"
#gwasR[gwasR$CHR=="4",2] <- "2A"
#gwasR[gwasR$CHR=="5",2] <- "2B"
#gwasR[gwasR$CHR=="6",2] <- "2D"
#gwasR[gwasR$CHR=="7",2] <- "3A"
#gwasR[gwasR$CHR=="8",2] <- "3B"
#gwasR[gwasR$CHR=="9",2] <- "3D"
#gwasR[gwasR$CHR=="10",2] <- "4A"
#gwasR[gwasR$CHR=="11",2] <- "4B"
#gwasR[gwasR$CHR=="12",2] <- "4D"
#gwasR[gwasR$CHR=="13",2] <- "5A"
#gwasR[gwasR$CHR=="14",2] <- "5B"
#gwasR[gwasR$CHR=="15",2] <- "5D"
#gwasR[gwasR$CHR=="16",2] <- "6A"
#gwasR[gwasR$CHR=="17",2] <- "6B"
#gwasR[gwasR$CHR=="18",2] <- "6D"
#gwasR[gwasR$CHR=="19",2] <- "7A"
#gwasR[gwasR$CHR=="20",2] <- "7B"
#gwasR[gwasR$CHR=="21",2] <- "7D"
#gwasR$CHR <- as.numeric(gwasR$CHR)
manhattan(gwasR, annotateTop = T, highlight = highsnp, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,25))
#manhattan(gwasR, main="Manhattan plot", ylim=c(0, 10), cex=0.6, cex.axis=0.9, col = c("blue","orange"), suggestiveline = F, genomewideline = F, chrlabs = c(1:21))
#qq(gwasR$P, main="Q-Q plot of GWAS p-value", xlim=c(0,7), ylim=c(0,12), pch=18, col = "blue4", cex=1.5, las=1)
#South_A
abline(h = 3.849, col="red", lwd=2, lty=2)
#South_B
abline(h = 3.73, col="red", lwd=2, lty=2)
#South_D
abline(h = 3.58, col="red", lwd=2, lty=2)
#North_A
abline(h = 3.76, col="red", lwd=2, lty=2)
#North_B
abline(h = 3.98, col="red", lwd=2, lty=2)
#North_D
abline(h = 3.61, col="red", lwd=2, lty=2)

data1 <- read.table("North_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#data2 <- data1[-which(data1$p == "NaN"),]
data1$SNP <- c(1:332899)
#gwasR <- data2[,c(2,3,4,7)]
gwasR <- data1[,c(7,1,2,5)]
#gwasR$P <- 10^(-gwasR$MeanY)
colnames(gwasR) <- c("SNP", "CHR", "BP","P")

highsnp <- c(99120,99232,99317,99324,99381,99530,99629,99658,99668,99693,99720,99752,99767,99808,99844,99847,99860,99888,99900,99917,99919,99922,99926,99941,99956)
manhattan(gwasR, annotateTop = T, highlight = highsnp, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,25))
