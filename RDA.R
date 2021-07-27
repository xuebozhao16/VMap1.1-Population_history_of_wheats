library(vegan)
library(ggplot2)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/climate_change")
taxa <- read.table("select_taxa.txt",header=T,stringsAsFactors = F)
taxa_EA <- taxa[which(taxa$Region=="EA"),1]
taxa_WA <- taxa[which(taxa$Region=="WA"),1]
taxa_SCA <- taxa[which(taxa$Region=="SCA"),1]
taxa_AF <- taxa[which(taxa$Region=="AF"),1]
taxa_EU <- taxa[which(taxa$Region=="EU"),1]

phylum <- read.delim('All_noMiss_0.05_2000.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
row.names(phylum) <- c(1:2000)
phylum <- data.frame(t(phylum))
env <- read.delim('select_bio.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env_all <- data.frame(env[,1:20])
env_temp <- env_all[,1:12]
env_prec <- env_all[,13:20]
#直接使用原始数据，不做转化。对于群落物种组成数据来讲（因为通常包含很多 0 值），不是很推荐
#rda_result <- rda(phylum~., env, scale = FALSE)
##tb-RDA

#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
phylum_hel <- decostand(phylum, method = 'hellinger')
#使用全部的环境数据
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb_temp <- rda(phylum_hel~., env_temp, scale = FALSE)
rda_tb_prec <- rda(phylum_hel~., env_prec, scale = FALSE)
#rda_tb.scaling1 <- summary(rda_tb, scaling = 2)
#rda_tb.scaling1
#若只关注局部环境数据，除了在原始表格中修改变量个数外，还可直接在 rda() 中指定
#rda_part <- rda(phylum~elevation+one+two+three+four+five+six+seven+eight+nine+ten+eleven+twelve+thirteen+fourteen+fifteen+sixteen+seventeen+eighteen+nineteen, data = env, scale = FALSE)

label <- read.table("select_taxa.txt", header=T, stringsAsFactors = F)
label$Region <- as.factor(label$Region)
label$cols  = label$Region 
label$cols  = gsub("EA","#97FFFF",label$cols)
label$cols  = gsub("EU","#FFE4E1",label$cols)
label$cols  = gsub("SCA","#FF8247",label$cols)
label$cols  = gsub("WA","#FF6A6A",label$cols)
label$cols  = gsub("AF","#D8BFD8",label$cols)
label$cols  = gsub("AM","#838B8B",label$cols)
cols <- label$cols

plot(rda_tb_all, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，双序图')
#points(rda_tb, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
points(rda_tb_all, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = alpha(cols, 0.7), cex = 1)

legend("topright", c("EA", "EU", "SCA","WA", "AF", "AM"),pch=19,col = c("#97FFFF","#FFE4E1","#FF8247","#FF6A6A","#D8BFD8","#838B8B"))
text(rda_tb_all, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

###现在是画选择信号
library(qqman)
library(tidyverse)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/climate_change/selection")
##SCA和南方的比较
data1 <- read.table("South_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#根据文件行数修改
data1$SNP <- c(1:330025)
gwasR <- data1[,c(7,1,2,5)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
manhattan(gwasR, annotateTop = T, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), 
          suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,30))
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

##SCA和北方的比较
data1 <- read.table("North_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#根据文件行数修改
data1$SNP <- c(1:332899)
gwasR <- data1[,c(7,1,2,5)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
manhattan(gwasR, annotateTop = T, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"),
          suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,30))


##GO分析
library(stringr)
library(ggplot2)
library(cowplot)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/climate_change/gene/window_region_gene")
GO_BP = read.table("Go_analyze_North.txt",header = T,sep="\t")
iris1 <- ggplot(data=GO_BP)+
  geom_bar(aes(x=Term,y=Count, fill=(FDR)), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(FDR),low="#7f518f", high = "#eedcf5") +
  ylab("Gene count") +
  expand_limits(y=c(0,0.8))+
  ylim(0,400)+
  theme(
    axis.text.x=element_text(color="black",size=rel(0.8)),
    axis.text.y=element_text(color="black", size=rel(1)),
    axis.title.x = element_text(color="black", size=rel(1)),
    axis.title.y = element_blank(),
    legend.text=element_text(color="black",size=rel(0.5)),
    legend.title = element_text(color="black",size=rel(0.7))
  ) +theme_bw() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
GO_BP = read.table("Go_analyze_South.txt",header = T,sep="\t")
iris2 <- ggplot(data=GO_BP)+
  geom_bar(aes(x=Term,y=Count, fill=(FDR)), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(FDR),low="#7f518f", high = "#eedcf5") +
  ylab("Gene count") +
  expand_limits(y=c(0,8))+
  ylim(0,400)+
  theme(
    axis.text.x=element_text(color="black",size=rel(0.8)),
    axis.text.y=element_text(color="black", size=rel(1)),
    axis.title.x = element_text(color="black", size=rel(1)),
    axis.title.y = element_blank(),
    legend.text=element_text(color="black",size=rel(0.5)),
    legend.title = element_text(color="black",size=rel(0.7))
    #legend.position=c(0,1),legend.justification=c(-1,0)
    #legend.position="top",
  )+theme_bw() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
plot_grid(iris1, iris2)



###下面的代码是查看已经克隆的基因在是不是有受到选择，以及受到选择的程度
library(heatmap.plus)
library(pheatmap)
library(RColorBrewer)
##A
display.brewer.all()
brewer.pal(9,'BuPu')
cols = c("#e9d5f0","#88419D")
data <- read.table("/Users/xuebozhao/Documents/LuLab/wheaintrogression/volcanofinder/GO/Gene_Go_A.txt",header=T,
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
















