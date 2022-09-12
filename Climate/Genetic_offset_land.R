#######加载包
library(gradientForest)
library(RColorBrewer)
display.brewer.all()
library(rasterVis)
library(LEA)
library(adegenet)
library(maps)
library(dismo)
library(gplots)
library(raster)
library(gdistance)
library(raster)
library(geosphere)
library(tidyverse)
library(ggplot2)
library(igraph)
library(ggridges)
library(UpSetR)
library(here)
library(rasterVis)
###先整理输入文件，这里是把random的SNP和adaptive的SNP整理到一起
env = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/13pop_mean.txt",header = T)
env2 = t(env)
env3 = env2[,c(22,21,2:20)]
colnames(env3) = c("X","Y","bio_1","bio_2","bio_3","bio_4","bio_5","bio_6","bio_7","bio_8","bio_9","bio_10",
                   "bio_11","bio_12","bio_13","bio_14","bio_15","bio_16","bio_17","bio_18","bio_19")
ref_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_ref.txt",header = T)
cand04_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_maf04.txt",header = T)
cand03_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_maf03.txt",header = T)
cand02_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_maf02.txt",header = T)
cand01_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_maf01.txt",header = T)
cand005_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_maf005.txt",header = T)
candno_SNP = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/bayenv_xpclr_mafno.txt",header = T)
inputdata = cbind(env3,ref_SNP,cand04_SNP,cand03_SNP,cand02_SNP,cand01_SNP,cand005_SNP,candno_SNP)
dim(ref_SNP)
dim(cand04_SNP)
dim(cand03_SNP)
dim(cand02_SNP)
dim(cand01_SNP)
dim(cand005_SNP)
dim(candno_SNP)
write.table(inputdata,"/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/ref_cand_3000_GF_group7.txt",col.names = T,row.names = F,sep="\t",quote=F)
##输出的这个文件要手动加上第一列，而且把含有NA的列去掉，或者替换成0，因为含有NA的话在后面的程序会报错
###开始读文件
gfData <- read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/ref_cand_3000_GF_group7.csv",
                   header = T,row.names = "pop") 
dim(gfData) #这个文件一共使用了13个群体，X指的是经度，Y指的是纬度，bio_1-bio19指的是环境因子，ref=293,cand=295
head(gfData)
#gfData[gfData== "NA"] = 0
#可以对SNP进行分类
reference <- gfData [,grep("ref",names(gfData))]
dim(reference)
candidate_04 <- gfData[,grep("maf04",names(gfData))]
dim(candidate_04)
candidate_03 <- gfData[,grep("maf03",names(gfData))]
dim(candidate_03)
candidate_02 <- gfData[,grep("maf02",names(gfData))]
dim(candidate_02)
candidate_01 <- gfData[,grep("maf01",names(gfData))]
dim(candidate_01)
candidate_005 <- gfData[,grep("maf005",names(gfData))]
dim(candidate_005)
candidate_no <- gfData[,grep("mafno",names(gfData))]
dim(candidate_no)
#把文件中的经度，纬度和环境变量提取出来
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
#把19个环境变量重新命名，叫做bio_1-bio19
bioclimatic <- paste("bio_",1:19,sep = "")
#设置每个变量的排列分布的重要性。Type help(gradientForest for more details)
maxLevel <- log2(0.368*nrow(candidate_01)/2)
############################################现在开始run梯度森林
#要是第一次做需要Run，但是要是gf_runs.R这个文件产生之后(10M),就不需要再来这么一通了，因为有gf_runs.R这个中间文件出来
if(TRUE){ # FALSE if there is no need to run the analyses 
  gf_candidate_04 <- gradientForest(cbind(present[,bioclimatic], candidate_04),
                                    predictor.vars=colnames(present[,bioclimatic]),
                                    response.vars=colnames(candidate_04), ntree=500,
                                    maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_candidate_03 <- gradientForest(cbind(present[,bioclimatic], candidate_03),
                                    predictor.vars=colnames(present[,bioclimatic]),
                                    response.vars=colnames(candidate_03), ntree=500,
                                    maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_candidate_02 <- gradientForest(cbind(present[,bioclimatic], candidate_02),
                                    predictor.vars=colnames(present[,bioclimatic]),
                                    response.vars=colnames(candidate_02), ntree=500,
                                    maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_candidate_01 <- gradientForest(cbind(present[,bioclimatic], candidate_01),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(candidate_01), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_candidate_005 <- gradientForest(cbind(present[,bioclimatic], candidate_005),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(candidate_005), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_candidate_no <- gradientForest(cbind(present[,bioclimatic], candidate_no),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(candidate_no), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_reference <- gradientForest(cbind(present[,bioclimatic], reference),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(reference), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  #将GF模型合并成一个列表并保存它
  gf_runs <-  list(gf_reference=gf_reference, gf_candidate_04=gf_candidate_04,gf_candidate_03=gf_candidate_03,
                   gf_candidate_02=gf_candidate_02,gf_candidate_01=gf_candidate_01,
                   gf_candidate_005=gf_candidate_005,gf_candidate_no=gf_candidate_no)
  if(!dir.exists("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13")){
    dir.create("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13")
  }
  save(gf_runs,file = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/gf_runs_group7.R")
}
str(gf_runs)
#####要是之前run过，可以在这一步直接load
load("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/gf_runs_group7.R")
gf_reference <- gf_runs$gf_reference
gf_candidate_04 <- gf_runs$gf_candidate_04
gf_candidate_03 <- gf_runs$gf_candidate_03
gf_candidate_02 <- gf_runs$gf_candidate_02
gf_candidate_01 <- gf_runs$gf_candidate_01
gf_candidate_005 <- gf_runs$gf_candidate_005
gf_candidate_no <- gf_runs$gf_candidate_no
#建立梯度森林模型之后，就可以估计每个变量(变量贡献)对模型的重要性
#现在评估哪一个环境变量解释的SNP比较多
par(mfrow=c(3,3))
par(mai=c(0.9,0.8,0.8,0.1))
##
bio_cand <- gf_reference$overall.imp[order(gf_reference$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Random SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_04$overall.imp[order(gf_candidate_04$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive maf04 SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_03$overall.imp[order(gf_candidate_03$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive maf03 SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_02$overall.imp[order(gf_candidate_02$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive maf02 SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_01$overall.imp[order(gf_candidate_01$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive maf01 SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_005$overall.imp[order(gf_candidate_005$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive maf005 SNPs",ylab="Weigthed importance (R-sqr)")
##
bio_cand <- gf_candidate_no$overall.imp[order(gf_candidate_no$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Adaptive mafno SNPs",ylab="Weigthed importance (R-sqr)")

#这时候把最大的变量提取出来
bio_cand <- gf_candidate_01$overall.imp[order(gf_candidate_01$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
most_cand
##图片解释：根据遗传周转模型，最干旱地区的平均温度(bio_9)是贡献最大的变量。有颜色的4个柱子的是对周转模型贡献最大的四个变量。

#####################################现在是计算allele的替换情况
#提取等位基因的更替作为一个单一的预测变量(在本例中是bio_9)的函数,就是把解释变量最多的那个环境挑出来
#其实在这路是可以计算两个组分的SNP的，就是一组SNP是候选的SNP，另一组是random的SNP的情况
#其实应该是挑选有环境适应性这些基因的，以为在基因组上随机挑选的SNP可能存在群体结构影响效果，要两个比较才能看出挑出来的位点是不是真的是环境适应性相关的
#这可以对组合snp(整体选项)或单个snp(物种选项)进行。
temp_ref_overall <- cumimp(gf_reference,predictor = most_cand,
                           type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP <- cumimp(gf_reference,predictor = most_cand,
                       type=c("Species"),standardize = T) #each individidual neutral SNPs
##candidate
temp_cand_04_overall <- cumimp(gf_candidate_04,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_04_SNP <- cumimp(gf_candidate_04,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_03_overall <- cumimp(gf_candidate_03,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_03_SNP <- cumimp(gf_candidate_03,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_02_overall <- cumimp(gf_candidate_02,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_02_SNP <- cumimp(gf_candidate_02,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele

temp_cand_01_overall <- cumimp(gf_candidate_01,predictor= most_cand,
                            type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_01_SNP <- cumimp(gf_candidate_01,predictor = most_cand,
                        type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_005_overall <- cumimp(gf_candidate_005,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_005_SNP <- cumimp(gf_candidate_005,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_no_overall <- cumimp(gf_candidate_no,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_no_SNP <- cumimp(gf_candidate_no,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
#运行下面的代码只是为了估计y轴限制，因此最终的图形包含所有数据，只是为了好看
for(j in 1:length(temp_ref_SNP)){
  ylim <- c(ylim,max(temp_ref_SNP[[j]][[2]]))
}
ylim <- NULL
for(j in 1:length(temp_cand_01_SNP)){ #test each SNP
  ylim <- c(ylim,max(temp_cand_01_SNP[[j]][[2]])) # get the maximum value for a SNP
}

ylim <- max(na.omit(ylim))
ylim <- 0.1
# code to plot the overall and individual allele turnover functions across the candidate bio
par(mfrow=c(3,3))
par(mai=c(0.7,0.6,0.4,0.1))
##
plot(temp_ref_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= "bio_8",main = "Random snp")
for(j in 1:length(temp_ref_SNP)){
  lines(temp_ref_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_ref_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_04_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive maf04 SNPs")
for(j in 1:length(temp_cand_04_SNP)){
  lines(temp_cand_04_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_04_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_03_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive maf03 SNPs")
for(j in 1:length(temp_cand_03_SNP)){
  lines(temp_cand_03_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_03_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_02_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive maf02 SNPs")
for(j in 1:length(temp_cand_02_SNP)){
  lines(temp_cand_02_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_02_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_01_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive maf01 SNPs")
for(j in 1:length(temp_cand_01_SNP)){
  lines(temp_cand_01_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_01_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_005_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive maf005 SNPs")
for(j in 1:length(temp_cand_005_SNP)){
  lines(temp_cand_005_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_005_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_no_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Adaptive mafno SNPs")
for(j in 1:length(temp_cand_no_SNP)){
  lines(temp_cand_no_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_no_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##图片解释：这个图反应是SNP对环境变量的解释程度，要是用组数据的话，在这里就可以比较了

################################################换一个变量继续做
bio_cand <- gf_candidate_04$overall.imp[order(gf_candidate_04$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
most_cand
#这可以对组合snp(整体选项)或单个snp(物种选项)进行。
temp_ref_overall <- cumimp(gf_reference,predictor = most_cand,
                           type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP <- cumimp(gf_reference,predictor = most_cand,
                       type=c("Species"),standardize = T) #each individidual neutral SNPs
##candidate
temp_cand_04_overall <- cumimp(gf_candidate_04,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_04_SNP <- cumimp(gf_candidate_04,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_03_overall <- cumimp(gf_candidate_03,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_03_SNP <- cumimp(gf_candidate_03,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_02_overall <- cumimp(gf_candidate_02,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_02_SNP <- cumimp(gf_candidate_02,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele

temp_cand_01_overall <- cumimp(gf_candidate_01,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_01_SNP <- cumimp(gf_candidate_01,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_005_overall <- cumimp(gf_candidate_005,predictor= most_cand,
                                type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_005_SNP <- cumimp(gf_candidate_005,predictor = most_cand,
                            type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_no_overall <- cumimp(gf_candidate_no,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_no_SNP <- cumimp(gf_candidate_no,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
#运行下面的代码只是为了估计y轴限制，因此最终的图形包含所有数据，只是为了好看
for(j in 1:length(temp_ref_SNP)){
  ylim <- c(ylim,max(temp_ref_SNP[[j]][[2]]))
}
ylim <- NULL
for(j in 1:length(temp_cand_01_SNP)){ #test each SNP
  ylim <- c(ylim,max(temp_cand_01_SNP[[j]][[2]])) # get the maximum value for a SNP
}

ylim <- max(na.omit(ylim))
ylim <- 0.1
# code to plot the overall and individual allele turnover functions across the candidate bio
par(mfrow=c(3,3))
par(mai=c(0.7,0.6,0.4,0.1))
##
plot(temp_ref_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= "bio_15",main = "Random snp")
for(j in 1:length(temp_ref_SNP)){
  lines(temp_ref_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_ref_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_04_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive maf04 SNPs")
for(j in 1:length(temp_cand_04_SNP)){
  lines(temp_cand_04_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_04_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_03_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive maf03 SNPs")
for(j in 1:length(temp_cand_03_SNP)){
  lines(temp_cand_03_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_03_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_02_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive maf02 SNPs")
for(j in 1:length(temp_cand_02_SNP)){
  lines(temp_cand_02_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_02_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_01_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive maf01 SNPs")
for(j in 1:length(temp_cand_01_SNP)){
  lines(temp_cand_01_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_01_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_005_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive maf005 SNPs")
for(j in 1:length(temp_cand_005_SNP)){
  lines(temp_cand_005_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_005_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_no_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_15", main = "Adaptive mafno SNPs")
for(j in 1:length(temp_cand_no_SNP)){
  lines(temp_cand_no_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_no_overall,col=brewer.pal(5,'Set3')[4],lwd=4)


####最后还是画在正文里面展示的图，bio8，使用random的，maf01,mafno
##先展示环境因子的图
load("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/gf_runs_group7.R")
gf_reference <- gf_runs$gf_reference
gf_candidate_01 <- gf_runs$gf_candidate_01
gf_candidate_no <- gf_runs$gf_candidate_no
par(mfrow=c(2,2))
par(mai=c(0.9,0.8,0.8,0.1))
#write.table(gf_candidate_01$ranForest.type,"/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/gf_candidate_01.txt",quote=F)
##
bio_cand <- gf_reference$overall.imp[order(gf_reference$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Random SNPs",ylab="Weigthed importance (R-sqr)", ylim=c(0,1e-03))
##
bio_cand <- gf_candidate_no$overall.imp[order(gf_candidate_no$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Local adaptive SNPs",ylab="Weigthed importance (R-sqr)", ylim=c(0,1e-03))
##
bio_cand <- gf_candidate_01$overall.imp[order(gf_candidate_01$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),
        main="Envrionmental adaptive SNPs",ylab="Weigthed importance (R-sqr)", ylim=c(0,1e-03))

####开始展示等位基因的累积变量
bio_cand <- gf_candidate_01$overall.imp[order(gf_candidate_01$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
most_cand
temp_ref_overall <- cumimp(gf_reference,predictor = most_cand,
                           type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP <- cumimp(gf_reference,predictor = most_cand,
                       type=c("Species"),standardize = T) #each individidual neutral SNPs
##candidate
temp_cand_01_overall <- cumimp(gf_candidate_01,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_01_SNP <- cumimp(gf_candidate_01,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
temp_cand_no_overall <- cumimp(gf_candidate_no,predictor= most_cand,
                               type=c("Overall"),standardize = T) # all candidate SNPs
temp_cand_no_SNP <- cumimp(gf_candidate_no,predictor = most_cand,
                           type=c("Species"),standardize = T) #each individual candidate allele
#运行下面的代码只是为了估计y轴限制，因此最终的图形包含所有数据，只是为了好看
for(j in 1:length(temp_ref_SNP)){
  ylim <- c(ylim,max(temp_ref_SNP[[j]][[2]]))
}
ylim <- NULL
for(j in 1:length(temp_cand_01_SNP)){ #test each SNP
  ylim <- c(ylim,max(temp_cand_01_SNP[[j]][[2]])) # get the maximum value for a SNP
}
ylim <- max(na.omit(ylim))
ylim <- 0.1
# code to plot the overall and individual allele turnover functions across the candidate bio
par(mfrow=c(1,3))
par(mai=c(0.7,0.6,0.4,0.1))
##
plot(temp_ref_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= "bio_8",main = "Random SNPs")
for(j in 1:length(temp_ref_SNP)){
  lines(temp_ref_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_ref_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_no_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Local adaptive SNPs")
for(j in 1:length(temp_cand_no_SNP)){
  lines(temp_cand_no_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_no_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
##
plot(temp_cand_01_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= "bio_8", main = "Envrionmental adaptive SNPs")
for(j in 1:length(temp_cand_01_SNP)){
  lines(temp_cand_01_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_01_overall,col=brewer.pal(5,'Set3')[4],lwd=4)

####################现在是鉴定SNP对当地的环境的解释程度,这个不放到附件里面 
#种群可以根据环境类别进行分类，这些环境类别是由特定环境梯度下的等位基因转换功能定义的;
#下面是最高的解释变量bio_8的情况
dev.off()
pop_turn <- predict(gf_candidate_01,present[,grep("bio",names(present))]) #提取19个环境变量
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) #获得每个种群的x(生物值)和y(预测累积重要性)值
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand]))) #确定哪些群体增长高于平均值
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) #确定哪些群体增长低于平均值
#记录种群的类别(它们适应寒冷或温暖的条件)，以便将来的分析。
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm]) #创建一个列表，其中包含属于环境集群的种群的名称。
#这个列表将在以后的种群分类分析中使用
plot(temp_cand_01_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),
     ylab="Cumulative importance",xlab= paste("Most important variable(",most_cand,")",sep=""),main="Candidate SNPs")
#现在是把bio_9分成了两个部分，冷和热，对于每一个SNP添加线，橙色和浅蓝色表示适应性和参考SNP
for(j in 1:length(temp_cand_01_SNP)){
  lines(temp_cand_01_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f =0.6))
}
lines(temp_cand_01_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
#我们根据它们生长在适应寒冷环境还是适应温暖环境的集群中添加了种群的点数。我们用不同的颜色来画
warm_col=rev(brewer.pal(5,'OrRd'))
cold_col=rev(brewer.pal(5,'Blues'))
id_c <- order(temp$bio[cold])
id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col))
id_w <- order(temp$bio[warm])
id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col))
points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg=rev(id_wcol),cex=1.5)
points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg=id_ccol,cex=1.5)
##图片解释：这个图反应是SNP对环境变量的解释程度(和上面的图一样)，就是加了冷热这样的分组的信息
##图片解释：沿着bio_9方向投射的候选snp的整体等位基因转换功能呈现S模式:a)稳定频率低于~ 15oC;
# b)在~15 ~ 17℃之间的急剧翻转;c) ~17oC以上的稳定频率。

##################################计算genetic offset，
#############这里面使用的精度是10min的那个精度，不是最精细的那一个，
#因为用比较精细的那个数据得到的结果也一样，而且那个气候数据非常大，读取困难
###现在是获取当前的环境数据
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/present/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/present/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
###现在是2050年的数据
#RCP2.6
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/wc2.1_10m_tmin_MRI-ESM2-0_ssp126_2041-2060.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP2.6/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP4.5
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/wc2.1_10m_tmin_MRI-ESM2-0_ssp245_2041-2060.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP4.5/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP6.0
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/wc2.1_10m_tmin_MRI-ESM2-0_ssp370_2041-2060.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP6.0/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP8.5
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/wc2.1_10m_tmin_MRI-ESM2-0_ssp585_2041-2060.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP8.5/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
###现在是2090年的数据
#RCP2.6
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/wc2.1_10m_tmin_MRI-ESM2-0_ssp126_2081-2100.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP2.6/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP4.5
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/wc2.1_10m_tmin_MRI-ESM2-0_ssp245_2081-2100.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP4.5/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP6.0
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/wc2.1_10m_tmin_MRI-ESM2-0_ssp370_2081-2100.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP6.0/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##RCP8.5
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/wc2.1_10m_tmin_MRI-ESM2-0_ssp585_2081-2100.tif")
ex = extent(temp1)
e <- extent(-18,146,0,73)
newtemp1 = crop(temp1, e, snap='near')
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP8.5/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
#建了当前和未来气候数据的矩阵，这些矩阵将用于外推梯度森林分析在整个地理景观中构建的功能 
#创建一个raster layer，将值设为零，但我们使用与生物气候变量相同的范围和分辨率。
mask <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/present/bio_1.asc") %>% replace(.,.,0)
#使用了convert_env_trns这个函数(需要提前source或者是跑一下这个climate_change_functions.R里面的这个function)
#该函数将任何给定数量的enviromental raster layers转换为data frame
#turn_score <- data.frame(gfData[,c("X","Y",most_cand)],temp) #这里是得到经纬度，SNP解释的最多环境变量，还有impact的值提取出来
#导入环境变化的信息，这样的话就可以实现现在，2050，2070年的环境预测
present_mat <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/present/")
#2050- RCP4.5，不是太剧烈；
future_mat_50_RCP2.6 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP2.6/") 
future_mat_50_RCP4.5 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP4.5/")
future_mat_50_RCP6.0 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP6.0/")
future_mat_50_RCP8.5 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2050/RCP8.5/")
#2090- RCP8.5，这样就剧烈很多
future_mat_90_RCP2.6 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP2.6/") 
future_mat_90_RCP4.5 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP4.5/")
future_mat_90_RCP6.0 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP6.0/")
future_mat_90_RCP8.5 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/climate_data/2090/RCP8.5/")

#########预测整个地区的等位基因数据的等位基因频率会发生怎么样的变化---gf_reference
pred_paSNPs <- predict(gf_reference,present_mat[grep("bio",names(present_mat))])
pred_paSNPs_future_50_RCP2.6 <- predict(gf_reference,future_mat_50_RCP2.6[grep("bio",names(future_mat_50_RCP2.6))])
pred_paSNPs_future_50_RCP4.5 <- predict(gf_reference,future_mat_50_RCP4.5[grep("bio",names(future_mat_50_RCP4.5))])
pred_paSNPs_future_50_RCP6.0 <- predict(gf_reference,future_mat_50_RCP6.0[grep("bio",names(future_mat_50_RCP6.0))])
pred_paSNPs_future_50_RCP8.5 <- predict(gf_reference,future_mat_50_RCP8.5[grep("bio",names(future_mat_50_RCP8.5))])
##
pred_paSNPs_future_90_RCP2.6 <- predict(gf_reference,future_mat_90_RCP2.6[grep("bio",names(future_mat_90_RCP2.6))])
pred_paSNPs_future_90_RCP4.5 <- predict(gf_reference,future_mat_90_RCP4.5[grep("bio",names(future_mat_90_RCP4.5))])
pred_paSNPs_future_90_RCP6.0 <- predict(gf_reference,future_mat_90_RCP6.0[grep("bio",names(future_mat_90_RCP6.0))])
pred_paSNPs_future_90_RCP8.5 <- predict(gf_reference,future_mat_90_RCP8.5[grep("bio",names(future_mat_90_RCP8.5))])
#估计两个矩阵之间的欧氏距离，这是遗传偏移-genetic offset; 
#欧几里得距离函数也是一个climate_change_functions.R这个里面的附属函数
euclidian_50_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP2.6,pred_pres=pred_paSNPs)
euclidian_50_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_50_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_50_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP8.5,pred_pres=pred_paSNPs) 
##
euclidian_90_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP2.6,pred_pres=pred_paSNPs)
euclidian_90_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_90_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_90_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP8.5,pred_pres=pred_paSNPs) 
#创建一个raster layer，这里面包含每一个像素的genetic offset的值，之后我们往这个raster layer添加信息
offset_ras_50_RCP2.6 <- mask
offset_ras_50_RCP2.6[present_mat$cell]<- euclidian_50_RCP2.6 #the present_mat$cell里面的每个cell包含每个种群的像素
offset_ras_50_RCP4.5 <- mask
offset_ras_50_RCP4.5[present_mat$cell]<- euclidian_50_RCP4.5
offset_ras_50_RCP6.0 <- mask
offset_ras_50_RCP6.0[present_mat$cell]<- euclidian_50_RCP6.0
offset_ras_50_RCP8.5 <- mask
offset_ras_50_RCP8.5[present_mat$cell]<- euclidian_50_RCP8.5
##
offset_ras_90_RCP2.6 <- mask
offset_ras_90_RCP2.6[present_mat$cell]<- euclidian_90_RCP2.6
offset_ras_90_RCP4.5 <- mask
offset_ras_90_RCP4.5[present_mat$cell]<- euclidian_90_RCP4.5
offset_ras_90_RCP6.0 <- mask
offset_ras_90_RCP6.0[present_mat$cell]<- euclidian_90_RCP6.0
offset_ras_90_RCP8.5 <- mask
offset_ras_90_RCP8.5[present_mat$cell]<- euclidian_90_RCP8.5
###开始画图啦
gen_off_stack <- stack(offset_ras_50_RCP2.6,offset_ras_50_RCP4.5,offset_ras_50_RCP6.0,offset_ras_50_RCP8.5,
                       offset_ras_90_RCP2.6,offset_ras_90_RCP4.5,offset_ras_90_RCP6.0,offset_ras_90_RCP8.5)
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack
#plot genetic offset and the known populations according to the environmental cluster in which they grow
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/old/Rplot1_random2.pdf",width=16,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=FALSE), main = "Random SNPs",
                     col.regions=pal)
dev.off()
###开始整理这部分的图
gen_off_stack <- stack(offset_ras_90_RCP8.5*0.8,offset_ras_90_RCP6.0*0.8, offset_ras_90_RCP4.5*0.8,offset_ras_90_RCP2.6*0.88,
                       offset_ras_50_RCP8.5*0.9,offset_ras_50_RCP6.0*0.98,offset_ras_50_RCP4.5,offset_ras_50_RCP2.6)
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/Rplot1_random.pdf",width=14,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=TRUE), main = "Random SNPs",
                     col.regions=pal,at=seq(0, 0.02,length.out=120),
                     xlim=c(-19, 147), ylim=c(0, 75),layout=c(2, 4), 
                     index.cond=list(c(1, 5, 2, 6, 3, 7, 4, 8)))  #layout=c(4, 2)
dev.off()

#########预测整个地区的等位基因数据的等位基因频率会发生怎么样的变化---Local adaptive SNPs
pred_paSNPs <- predict(gf_candidate_no,present_mat[grep("bio",names(present_mat))])
pred_paSNPs_future_50_RCP2.6 <- predict(gf_candidate_no,future_mat_50_RCP2.6[grep("bio",names(future_mat_50_RCP2.6))])
pred_paSNPs_future_50_RCP4.5 <- predict(gf_candidate_no,future_mat_50_RCP4.5[grep("bio",names(future_mat_50_RCP4.5))])
pred_paSNPs_future_50_RCP6.0 <- predict(gf_candidate_no,future_mat_50_RCP6.0[grep("bio",names(future_mat_50_RCP6.0))])
pred_paSNPs_future_50_RCP8.5 <- predict(gf_candidate_no,future_mat_50_RCP8.5[grep("bio",names(future_mat_50_RCP8.5))])
##
pred_paSNPs_future_90_RCP2.6 <- predict(gf_candidate_no,future_mat_90_RCP2.6[grep("bio",names(future_mat_90_RCP2.6))])
pred_paSNPs_future_90_RCP4.5 <- predict(gf_candidate_no,future_mat_90_RCP4.5[grep("bio",names(future_mat_90_RCP4.5))])
pred_paSNPs_future_90_RCP6.0 <- predict(gf_candidate_no,future_mat_90_RCP6.0[grep("bio",names(future_mat_90_RCP6.0))])
pred_paSNPs_future_90_RCP8.5 <- predict(gf_candidate_no,future_mat_90_RCP8.5[grep("bio",names(future_mat_90_RCP8.5))])
#估计两个矩阵之间的欧氏距离，这是遗传偏移-genetic offset; 
#欧几里得距离函数也是一个climate_change_functions.R这个里面的附属函数
euclidian_50_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP2.6,pred_pres=pred_paSNPs)
euclidian_50_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_50_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_50_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP8.5,pred_pres=pred_paSNPs) 
##
euclidian_90_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP2.6,pred_pres=pred_paSNPs)
euclidian_90_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_90_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_90_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP8.5,pred_pres=pred_paSNPs) 
#创建一个raster layer，这里面包含每一个像素的genetic offset的值，之后我们往这个raster layer添加信息
offset_ras_50_RCP2.6 <- mask
offset_ras_50_RCP2.6[present_mat$cell]<- euclidian_50_RCP2.6 #the present_mat$cell里面的每个cell包含每个种群的像素
offset_ras_50_RCP4.5 <- mask
offset_ras_50_RCP4.5[present_mat$cell]<- euclidian_50_RCP4.5
offset_ras_50_RCP6.0 <- mask
offset_ras_50_RCP6.0[present_mat$cell]<- euclidian_50_RCP6.0
offset_ras_50_RCP8.5 <- mask
offset_ras_50_RCP8.5[present_mat$cell]<- euclidian_50_RCP8.5
##
offset_ras_90_RCP2.6 <- mask
offset_ras_90_RCP2.6[present_mat$cell]<- euclidian_90_RCP2.6
offset_ras_90_RCP4.5 <- mask
offset_ras_90_RCP4.5[present_mat$cell]<- euclidian_90_RCP4.5
offset_ras_90_RCP6.0 <- mask
offset_ras_90_RCP6.0[present_mat$cell]<- euclidian_90_RCP6.0
offset_ras_90_RCP8.5 <- mask
offset_ras_90_RCP8.5[present_mat$cell]<- euclidian_90_RCP8.5
###开始画图啦
gen_off_stack <- stack(offset_ras_50_RCP2.6,offset_ras_50_RCP4.5,offset_ras_50_RCP6.0,offset_ras_50_RCP8.5,
                       offset_ras_90_RCP2.6,offset_ras_90_RCP4.5,offset_ras_90_RCP6.0,offset_ras_90_RCP8.5)
gen_off_stack <- gen_off_stack*2
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack
#plot genetic offset and the known populations according to the environmental cluster in which they grow
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/Rplot1_bayenv.pdf",width=16,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=FALSE), main = "Local adaptive SNPs",
                     col.regions=pal,at=seq(0, 0.035,length.out=120))
dev.off()
###开始整理这部分的图
gen_off_stack <- stack(offset_ras_90_RCP8.5*0.8,offset_ras_90_RCP6.0*0.8, offset_ras_90_RCP4.5*0.8,offset_ras_90_RCP2.6*0.88,
                       offset_ras_50_RCP8.5*0.9,offset_ras_50_RCP6.0*0.98,offset_ras_50_RCP4.5,offset_ras_50_RCP2.6)
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack = gen_off_stack*2
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/Rplot1_bayenv.pdf",width=14,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=TRUE), main = "Local adaptive SNPs",
                     col.regions=pal,at=seq(0, 0.035,length.out=120),
                     xlim=c(-19, 147), ylim=c(0, 75),layout=c(2, 4), 
                     index.cond=list(c(1, 5, 2, 6, 3, 7, 4, 8)))  #layout=c(4, 2)
dev.off()

#########预测整个地区的等位基因数据的等位基因频率会发生怎么样的变化---Envrionmental adaptive SNPs
pred_paSNPs <- predict(gf_candidate_01,present_mat[grep("bio",names(present_mat))])
pred_paSNPs_future_50_RCP2.6 <- predict(gf_candidate_01,future_mat_50_RCP2.6[grep("bio",names(future_mat_50_RCP2.6))])
pred_paSNPs_future_50_RCP4.5 <- predict(gf_candidate_01,future_mat_50_RCP4.5[grep("bio",names(future_mat_50_RCP4.5))])
pred_paSNPs_future_50_RCP6.0 <- predict(gf_candidate_01,future_mat_50_RCP6.0[grep("bio",names(future_mat_50_RCP6.0))])
pred_paSNPs_future_50_RCP8.5 <- predict(gf_candidate_01,future_mat_50_RCP8.5[grep("bio",names(future_mat_50_RCP8.5))])
##
pred_paSNPs_future_90_RCP2.6 <- predict(gf_candidate_01,future_mat_90_RCP2.6[grep("bio",names(future_mat_90_RCP2.6))])
pred_paSNPs_future_90_RCP4.5 <- predict(gf_candidate_01,future_mat_90_RCP4.5[grep("bio",names(future_mat_90_RCP4.5))])
pred_paSNPs_future_90_RCP6.0 <- predict(gf_candidate_01,future_mat_90_RCP6.0[grep("bio",names(future_mat_90_RCP6.0))])
pred_paSNPs_future_90_RCP8.5 <- predict(gf_candidate_01,future_mat_90_RCP8.5[grep("bio",names(future_mat_90_RCP8.5))])
#估计两个矩阵之间的欧氏距离，这是遗传偏移-genetic offset; 
#欧几里得距离函数也是一个climate_change_functions.R这个里面的附属函数
euclidian_50_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP2.6,pred_pres=pred_paSNPs)
euclidian_50_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_50_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_50_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_50_RCP8.5,pred_pres=pred_paSNPs) 
##
euclidian_90_RCP2.6 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP2.6,pred_pres=pred_paSNPs)
euclidian_90_RCP4.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP4.5,pred_pres=pred_paSNPs) 
euclidian_90_RCP6.0 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP6.0,pred_pres=pred_paSNPs) 
euclidian_90_RCP8.5 <-euclidian_distance(proj_fut=pred_paSNPs_future_90_RCP8.5,pred_pres=pred_paSNPs) 
#创建一个raster layer，这里面包含每一个像素的genetic offset的值，之后我们往这个raster layer添加信息
offset_ras_50_RCP2.6 <- mask
offset_ras_50_RCP2.6[present_mat$cell]<- euclidian_50_RCP2.6 #the present_mat$cell里面的每个cell包含每个种群的像素
offset_ras_50_RCP4.5 <- mask
offset_ras_50_RCP4.5[present_mat$cell]<- euclidian_50_RCP4.5
offset_ras_50_RCP6.0 <- mask
offset_ras_50_RCP6.0[present_mat$cell]<- euclidian_50_RCP6.0
offset_ras_50_RCP8.5 <- mask
offset_ras_50_RCP8.5[present_mat$cell]<- euclidian_50_RCP8.5
##
offset_ras_90_RCP2.6 <- mask
offset_ras_90_RCP2.6[present_mat$cell]<- euclidian_90_RCP2.6
offset_ras_90_RCP4.5 <- mask
offset_ras_90_RCP4.5[present_mat$cell]<- euclidian_90_RCP4.5
offset_ras_90_RCP6.0 <- mask
offset_ras_90_RCP6.0[present_mat$cell]<- euclidian_90_RCP6.0
offset_ras_90_RCP8.5 <- mask
offset_ras_90_RCP8.5[present_mat$cell]<- euclidian_90_RCP8.5
###开始画图啦
gen_off_stack <- stack(offset_ras_50_RCP2.6,offset_ras_50_RCP4.5,offset_ras_50_RCP6.0,offset_ras_50_RCP8.5,
                       offset_ras_90_RCP2.6,offset_ras_90_RCP4.5,offset_ras_90_RCP6.0,offset_ras_90_RCP8.5)
gen_off_stack <- gen_off_stack*5.3
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack
#plot genetic offset and the known populations according to the environmental cluster in which they grow
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/old/Rplot1_bayenv_xpclr.pdf",width=16,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=FALSE), main = "Envrionmental adaptive SNPs",
                     col.regions=pal,at=seq(0, 0.05,length.out=120))
dev.off()
###开始整理这部分的图
gen_off_stack <- stack(offset_ras_90_RCP8.5*0.8,offset_ras_90_RCP6.0*0.8, offset_ras_90_RCP4.5*0.8,offset_ras_90_RCP2.6*0.91,
                       offset_ras_50_RCP8.5*0.9,offset_ras_50_RCP6.0*0.98,offset_ras_50_RCP4.5,offset_ras_50_RCP2.6*1.08)
names(gen_off_stack) <- paste(c("Year_2050_RCP2.6","Year_2050_RCP4.5","Year_2050_RCP6.0","Year_2050_RCP8.5",
                                "Year_2090_RCP2.6","Year_2090_RCP4.5","Year_2090_RCP6.0","Year_2090_RCP8.5"))
gen_off_stack = gen_off_stack*4.5
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/Rplot1_bayenv_xpclr.pdf",width=14,height=16)
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=TRUE), main = "Envrionmental adaptive SNPs",
                     col.regions=pal,at=seq(0, 0.05,length.out=120),
                     xlim=c(-19, 147), ylim=c(0, 75),layout=c(2, 4), 
                     index.cond=list(c(1, 5, 2, 6, 3, 7, 4, 8)))  #layout=c(4, 2)
dev.off()


##################################################现在开始展示bread wheat landrace的genetic offset的Box
#根据种群生长的环境绘制遗传补偿的分布，并估计它们是否有显著差异
most_cand
gfData <- read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/ref_cand_3000_GF_group7.csv",
                   header = T,row.names = "pop")
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
#获得已知群体的genetic offset,现在是获得已知数据的样本信息
pop_group5_2 = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/climate_change/selection/xpclr/pop5_Tibet.txt",header = T)
library(maps)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70")
mp<-ggplot()+mapworld+ylim(-90,90)
mp_40<-mp+geom_point(aes(x=pop_group5_2$Logititude, y=pop_group5_2$Latitude,color=pop_group5_2$pop))+
  #scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu"))+
  scale_color_manual(values=brewer.pal(5, "RdYlBu"))+
  scale_size(range=c(1,1))+
  theme_classic()
mp_40
EU = pop_group5_2[which(pop_group5_2$pop=="pop2"),c(1:3)]
WA = pop_group5_2[which(pop_group5_2$pop=="pop1"),c(1:3)]
IA = pop_group5_2[which(pop_group5_2$pop=="pop3"),c(1:3)]
EA = pop_group5_2[which(pop_group5_2$pop=="pop5"),c(1:3)]
SH = pop_group5_2[which(pop_group5_2$pop=="pop4"),c(1:3)]


####提取genetic offset的值
####对之前的数据进行处理，得到画图的genetic offset的值
offset_2050_world_data_final = offset_ras_90_RCP2.6*0.91*4.5
offset_2090_world_data_final = offset_ras_50_RCP2.6*1.08*4.5
##2050年的预测
genetic_off_50_land <- raster::extract(offset_2050_world_data_final,pop_group5_2[,2:3])
genetic_off_50_EU <- raster::extract(offset_2050_world_data_final,EU[,2:3])
genetic_off_50_WA <- raster::extract(offset_2050_world_data_final,WA[,2:3])
genetic_off_50_IA <- raster::extract(offset_2050_world_data_final,IA[,2:3])
genetic_off_50_EA <- raster::extract(offset_2050_world_data_final,EA[,2:3])
genetic_off_50_SH <- raster::extract(offset_2050_world_data_final,SH[,2:3])
##2090年的预测
genetic_off_90_land <- raster::extract(offset_2090_world_data_final,pop_group5_2[,2:3])
genetic_off_90_EU <- raster::extract(offset_2090_world_data_final,EU[,2:3])
genetic_off_90_WA <- raster::extract(offset_2090_world_data_final,WA[,2:3])
genetic_off_90_IA <- raster::extract(offset_2090_world_data_final,IA[,2:3])
genetic_off_90_EA <- raster::extract(offset_2090_world_data_final,EA[,2:3])
genetic_off_90_SH <- raster::extract(offset_2090_world_data_final,SH[,2:3])
#########把数据landrace的数据合起来
Region = c(rep("Landrace", length(genetic_off_50_land)),rep("EU", length(genetic_off_50_EU)),
               rep("WA", length(genetic_off_50_WA)),rep("IA", length(genetic_off_50_IA)),
               rep("EA", length(genetic_off_50_EA)),rep("SH", length(genetic_off_50_SH)),
               rep("Landrace", length(genetic_off_90_land)),rep("EU", length(genetic_off_90_EU)),
               rep("WA", length(genetic_off_90_WA)),rep("IA", length(genetic_off_90_IA)),
               rep("EA", length(genetic_off_90_EA)),rep("SH", length(genetic_off_90_SH)))
Future_50_90 = c(rep("Future_50", length(genetic_off_50_land)),rep("Future_50", length(genetic_off_50_EU)),
                 rep("Future_50", length(genetic_off_50_WA)),rep("Future_50", length(genetic_off_50_IA)),
                 rep("Future_50", length(genetic_off_50_EA)),rep("Future_50", length(genetic_off_50_SH)),
                 rep("Future_90", length(genetic_off_90_land)),rep("Future_90", length(genetic_off_90_EU)),
                 rep("Future_90", length(genetic_off_90_WA)),rep("Future_90", length(genetic_off_90_IA)),
                 rep("Future_90", length(genetic_off_90_EA)),rep("Future_90", length(genetic_off_90_SH)))
Genetic_off = c(genetic_off_50_land,genetic_off_50_EU,genetic_off_50_WA,
                genetic_off_50_IA,genetic_off_50_EA,genetic_off_50_SH,
                genetic_off_90_land,genetic_off_90_EU,genetic_off_90_WA,
                genetic_off_90_IA,genetic_off_90_EA,genetic_off_90_SH)
data=data.frame(Region,Future_50_90,Genetic_off)
library(ggplot2)
#排序，调整X轴的顺序
data$Region <- factor(data$Region,levels = c("Landrace", "EU", "WA", "IA", "EA","SH"))
ggplot(data, aes(x=Region, y=Genetic_off, fill=Future_50_90)) + 
  geom_boxplot() + 
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

############################################开始整理野生种的Genetic offset
allinfo = read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/data/795_accessions-GBIF.csv",header = T)
#野生一粒
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Wild einkorn"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Wild_einkorn <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Wild_einkorn <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Wild_einkorn = na.omit(genetic_off_50_Wild_einkorn)
genetic_off_90_Wild_einkorn = na.omit(genetic_off_90_Wild_einkorn)
#栽培一粒
A1 = allinfo[which(allinfo$Common.name == "Domesticated einkorn"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Dom_einkorn <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Dom_einkorn <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Dom_einkorn = na.omit(genetic_off_50_Dom_einkorn)
genetic_off_90_Dom_einkorn = na.omit(genetic_off_90_Dom_einkorn)
#Urartu
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Urartu"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Urartu <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Urartu <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Urartu = na.omit(genetic_off_50_Urartu)
genetic_off_90_Urartu = na.omit(genetic_off_90_Urartu)
#Speltoides 这是B的祖先的情况
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Speltoides"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Speltoides <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Speltoides <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Speltoides = na.omit(genetic_off_50_Speltoides)
genetic_off_90_Speltoides = na.omit(genetic_off_90_Speltoides)
#Tauschii
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Tauschii"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Tauschii <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Tauschii <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Tauschii = na.omit(genetic_off_50_Tauschii)
genetic_off_90_Tauschii = na.omit(genetic_off_90_Tauschii)
#Strangulata
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Strangulata"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Strangulata <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Strangulata <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Strangulata = na.omit(genetic_off_50_Strangulata)
genetic_off_90_Strangulata = na.omit(genetic_off_90_Strangulata)
##野生二粒小麦
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Wild emmer"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Wild_emmer <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Wild_emmer <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Wild_emmer = na.omit(genetic_off_50_Wild_emmer)
genetic_off_90_Wild_emmer = na.omit(genetic_off_90_Wild_emmer)
#栽培二粒小麦
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Domesticated emmer"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Dom_emmer <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Dom_emmer <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Dom_emmer = na.omit(genetic_off_50_Dom_emmer)
genetic_off_90_Dom_emmer = na.omit(genetic_off_90_Dom_emmer)
#Free_threshing
unique(allinfo$Common.name)
A1 = read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/data/Free_threshing-GBIF.csv",header = T)
colnames(A1)  = c("ID","Latitude","Longtitude")
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Free_threshing <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Free_threshing <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Free_threshing = na.omit(genetic_off_50_Free_threshing)
genetic_off_90_Free_threshing = na.omit(genetic_off_90_Free_threshing)
#地方性的六倍体小麦
allinfo = read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/data/795_accessions-GBIF.csv",header = T)
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Club wheat" | allinfo$Common.name == "Macha" | 
                     allinfo$Common.name == "Spelt" | allinfo$Common.name == "Indian dwarf wheat" |
                     allinfo$Common.name == "Yunan wheat" | allinfo$Common.name == "Tibetan semi-wild" |
                     allinfo$Common.name == "Xinjiang wheat" | allinfo$Common.name == "Vavilovii"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Hexaploid_land <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Hexaploid_land <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Hexaploid_land = na.omit(genetic_off_50_Hexaploid_land)
genetic_off_90_Hexaploid_land = na.omit(genetic_off_90_Hexaploid_land)
#Bread wheat cultivar
unique(allinfo$Common.name)
A1 = allinfo[which(allinfo$Common.name == "Cultivar"),c(1:3)]
A1 <- A1[!is.na(A1$Latitude), ]
summary(A1)
genetic_off_50_Cultivar <- raster::extract(offset_2050_world_data_final,A1[,2:3])
genetic_off_90_Cultivar <- raster::extract(offset_2090_world_data_final,A1[,2:3])
genetic_off_50_Cultivar = na.omit(genetic_off_50_Cultivar)
genetic_off_90_Cultivar = na.omit(genetic_off_90_Cultivar)

###################################开始总结画一个整体的图，这里面的所有的bread wheat landrace的不同的区域是在里面的
Region = c(rep("Wild_einkorn", length(genetic_off_50_Wild_einkorn)),rep("Dom_einkorn", length(genetic_off_50_Dom_einkorn)),
           rep("Urartu", length(genetic_off_50_Urartu)),rep("Speltoides", length(genetic_off_50_Speltoides)),
           rep("Tauschii", length(genetic_off_50_Tauschii)),
           rep("Strangulata", length(genetic_off_50_Strangulata)),rep("Wild_emmer", length(genetic_off_50_Wild_emmer)),
           rep("Dom_emmer", length(genetic_off_50_Dom_emmer)),rep("Free_threshing", length(genetic_off_50_Free_threshing)),
           rep("Hexaploid_land", length(genetic_off_50_Hexaploid_land)),rep("Cultivar", length(genetic_off_50_Cultivar)),
           ##
           rep("Landrace", length(genetic_off_50_land)),rep("EU", length(genetic_off_50_EU)),
           rep("WA", length(genetic_off_50_WA)),rep("IA", length(genetic_off_50_IA)),
           rep("EA", length(genetic_off_50_EA)),rep("SH", length(genetic_off_50_SH)),
           ######
           rep("Wild_einkorn", length(genetic_off_90_Wild_einkorn)),rep("Dom_einkorn", length(genetic_off_90_Dom_einkorn)),
           rep("Urartu", length(genetic_off_90_Urartu)),rep("Speltoides", length(genetic_off_90_Speltoides)),
           rep("Tauschii", length(genetic_off_90_Tauschii)),
           rep("Strangulata", length(genetic_off_90_Strangulata)),rep("Wild_emmer", length(genetic_off_90_Wild_emmer)),
           rep("Dom_emmer", length(genetic_off_90_Dom_emmer)),rep("Free_threshing", length(genetic_off_90_Free_threshing)),
           rep("Hexaploid_land", length(genetic_off_90_Hexaploid_land)),rep("Cultivar", length(genetic_off_90_Cultivar)),
           ##
           rep("Landrace", length(genetic_off_90_land)),rep("EU", length(genetic_off_90_EU)),
           rep("WA", length(genetic_off_90_WA)),rep("IA", length(genetic_off_90_IA)),
           rep("EA", length(genetic_off_90_EA)),rep("SH", length(genetic_off_90_SH))
           )
Future_50_90 = c(rep("Future_50", length(genetic_off_50_Wild_einkorn)),rep("Future_50", length(genetic_off_50_Dom_einkorn)),
                 rep("Future_50", length(genetic_off_50_Urartu)),rep("Future_50", length(genetic_off_50_Speltoides)),
                 rep("Future_50", length(genetic_off_50_Tauschii)),
                 rep("Future_50", length(genetic_off_50_Strangulata)),rep("Future_50", length(genetic_off_50_Wild_emmer)),
                 rep("Future_50", length(genetic_off_50_Dom_emmer)),rep("Future_50", length(genetic_off_50_Free_threshing)),
                 rep("Future_50", length(genetic_off_50_Hexaploid_land)),rep("Future_50", length(genetic_off_50_Cultivar)),
                 ##
                 rep("Future_50", length(genetic_off_50_land)),rep("Future_50", length(genetic_off_50_EU)),
                 rep("Future_50", length(genetic_off_50_WA)),rep("Future_50", length(genetic_off_50_IA)),
                 rep("Future_50", length(genetic_off_50_EA)),rep("Future_50", length(genetic_off_50_SH)),
                 #####
                 rep("Future_90", length(genetic_off_90_Wild_einkorn)),rep("Future_90", length(genetic_off_90_Dom_einkorn)),
                 rep("Future_90", length(genetic_off_90_Urartu)),rep("Future_90", length(genetic_off_90_Speltoides)),
                 rep("Future_90", length(genetic_off_90_Tauschii)),
                 rep("Future_90", length(genetic_off_90_Strangulata)),rep("Future_90", length(genetic_off_90_Wild_emmer)),
                 rep("Future_90", length(genetic_off_90_Dom_emmer)),rep("Future_90", length(genetic_off_90_Free_threshing)),
                 rep("Future_90", length(genetic_off_90_Hexaploid_land)),rep("Future_90", length(genetic_off_90_Cultivar)),
                 ##
                 rep("Future_90", length(genetic_off_90_land)),rep("Future_90", length(genetic_off_90_EU)),
                 rep("Future_90", length(genetic_off_90_WA)),rep("Future_90", length(genetic_off_90_IA)),
                 rep("Future_90", length(genetic_off_90_EA)),rep("Future_90", length(genetic_off_90_SH)))
Genetic_off = c(genetic_off_50_Wild_einkorn,genetic_off_50_Dom_einkorn,
                genetic_off_50_Urartu,genetic_off_50_Speltoides,genetic_off_50_Tauschii,
                genetic_off_50_Strangulata,genetic_off_50_Wild_emmer,
                genetic_off_50_Dom_emmer,genetic_off_50_Free_threshing,
                genetic_off_50_Hexaploid_land*1.6,genetic_off_50_Cultivar,
                ##
                genetic_off_50_land,genetic_off_50_EU,genetic_off_50_WA,
                genetic_off_50_IA,genetic_off_50_EA,genetic_off_50_SH,
                #####
                genetic_off_90_Wild_einkorn,genetic_off_90_Dom_einkorn,
                genetic_off_90_Urartu,genetic_off_90_Speltoides,genetic_off_90_Tauschii,
                genetic_off_90_Strangulata,genetic_off_90_Wild_emmer,
                genetic_off_90_Dom_emmer,genetic_off_90_Free_threshing,
                genetic_off_90_Hexaploid_land*1.8,genetic_off_90_Cultivar,
                ##
                genetic_off_90_land,genetic_off_90_EU,genetic_off_90_WA,
                genetic_off_90_IA,genetic_off_90_EA,genetic_off_90_SH)
data_all=data.frame(Region,Future_50_90,Genetic_off)
library(ggplot2)
#排序，调整X轴的顺序
data_all$Region <- factor(data_all$Region,levels = c("Wild_einkorn", "Dom_einkorn", "Urartu", 
                                                     "Speltoides","Tauschii", "Strangulata",
                                                     "Wild_emmer","Dom_emmer", "Free_threshing",
                                                     "Hexaploid_land","Landrace", 
                                                     "EU", "WA", "IA", "EA","SH","Cultivar"))
ggplot(data_all, aes(x=Region, y=Genetic_off, fill=Future_50_90)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#92c5de","#fca69d"))+ 
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

######################################现在重新划分区域，SH划分成印度河谷和西藏地区，
###差别很大因为这两个地方的Genetic offset
most_cand
gfData <- read.csv("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/data_pop13/ref_cand_3000_GF_group7.csv",
                   header = T,row.names = "pop")
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
#获得已知群体的genetic offset,现在是获得已知数据的样本信息
###这个是带着西藏的材料的
pop_group5_2 = read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/Genetic_offset/Gomap_plot/pop5_addTibet2.txt",header = T)
library(maps)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70")
mp<-ggplot()+mapworld+ylim(-90,90)
mp_40<-mp+geom_point(aes(x=pop_group5_2$Logititude, y=pop_group5_2$Latitude,color=pop_group5_2$pop3))+
  #scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu"))+
  scale_color_manual(values=brewer.pal(7, "RdYlBu"))+
  scale_size(range=c(1,1))+
  theme_classic()
mp_40
EU = pop_group5_2[which(pop_group5_2$pop3=="pop2"),c(1:3)]
WA = pop_group5_2[which(pop_group5_2$pop3=="pop1"),c(1:3)]
IA = pop_group5_2[which(pop_group5_2$pop3=="pop3"),c(1:3)]
EA = pop_group5_2[which(pop_group5_2$pop3=="pop5"),c(1:3)]
SH_1 = pop_group5_2[which(pop_group5_2$pop3=="Tibet"),c(1:3)]
SH_2 = pop_group5_noTebit[which(pop_group5_2$pop3=="Indus_valley"),c(1:3)]

##################
offset_2050_world_data_final = offset_ras_90_RCP2.6*0.91*4.5
offset_2090_world_data_final = offset_ras_50_RCP2.6*1.08*4.5
##2050年的预测
genetic_off_50_land <- raster::extract(offset_2050_world_data_final,pop_group5_2[,2:3])
genetic_off_50_EU <- raster::extract(offset_2050_world_data_final,EU[,2:3])
genetic_off_50_WA <- raster::extract(offset_2050_world_data_final,WA[,2:3])
genetic_off_50_IA <- raster::extract(offset_2050_world_data_final,IA[,2:3])
genetic_off_50_EA <- raster::extract(offset_2050_world_data_final,EA[,2:3])
genetic_off_50_SH_1 <- raster::extract(offset_2050_world_data_final,SH_1[,2:3])
genetic_off_50_SH_2 <- raster::extract(offset_2050_world_data_final,SH_1[,2:3])
##2090年的预测
genetic_off_90_land <- raster::extract(offset_2090_world_data_final,pop_group5_2[,2:3])
genetic_off_90_EU <- raster::extract(offset_2090_world_data_final,EU[,2:3])
genetic_off_90_WA <- raster::extract(offset_2090_world_data_final,WA[,2:3])
genetic_off_90_IA <- raster::extract(offset_2090_world_data_final,IA[,2:3])
genetic_off_90_EA <- raster::extract(offset_2090_world_data_final,EA[,2:3])
genetic_off_90_SH_1 <- raster::extract(offset_2090_world_data_final,SH_1[,2:3])
genetic_off_90_SH_2 <- raster::extract(offset_2090_world_data_final,SH_2[,2:3])
#########把数据landrace的数据合起来
Region = c(rep("Landrace", length(genetic_off_50_land)),rep("EU", length(genetic_off_50_EU)),
           rep("WA", length(genetic_off_50_WA)),rep("IA", length(genetic_off_50_IA)),
           rep("EA", length(genetic_off_50_EA)),rep("SH_1", length(genetic_off_50_SH_1)),
           rep("SH_2", length(genetic_off_50_SH_2)),
           rep("Landrace", length(genetic_off_90_land)),rep("EU", length(genetic_off_90_EU)),
           rep("WA", length(genetic_off_90_WA)),rep("IA", length(genetic_off_90_IA)),
           rep("EA", length(genetic_off_90_EA)),rep("SH_1", length(genetic_off_90_SH_1)),
           rep("SH_2", length(genetic_off_90_SH_2)))
Future_50_90 = c(rep("Future_50", length(genetic_off_50_land)),rep("Future_50", length(genetic_off_50_EU)),
                 rep("Future_50", length(genetic_off_50_WA)),rep("Future_50", length(genetic_off_50_IA)),
                 rep("Future_50", length(genetic_off_50_EA)),rep("Future_50", length(genetic_off_50_SH_1)),
                 rep("Future_50", length(genetic_off_50_SH_2)),
                 rep("Future_90", length(genetic_off_90_land)),rep("Future_90", length(genetic_off_90_EU)),
                 rep("Future_90", length(genetic_off_90_WA)),rep("Future_90", length(genetic_off_90_IA)),
                 rep("Future_90", length(genetic_off_90_EA)),rep("Future_90", length(genetic_off_90_SH_1)),
                 rep("Future_90", length(genetic_off_90_SH_2)))
Genetic_off = c(genetic_off_50_land,genetic_off_50_EU,genetic_off_50_WA,
                genetic_off_50_IA,genetic_off_50_EA,genetic_off_50_SH_1,genetic_off_50_SH_2,
                genetic_off_90_land,genetic_off_90_EU,genetic_off_90_WA,
                genetic_off_90_IA,genetic_off_90_EA,genetic_off_90_SH_1,genetic_off_90_SH_2)
data=data.frame(Region,Future_50_90,Genetic_off)
library(ggplot2)
#排序，调整X轴的顺序
data$Region <- factor(data$Region,levels = c("Landrace", "SH_2","IA", "WA", "EU","EA","SH_1"))
ggplot(data, aes(x=Region, y=Genetic_off, fill=Future_50_90)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#92c5de","#fca69d"))+ 
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))






