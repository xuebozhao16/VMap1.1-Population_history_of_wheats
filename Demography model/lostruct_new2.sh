##local一共有3步
#1.Compute the local PCA coordinates - done with eigen_windows(). 这个方法是为了得到一个numeric matrix，每一列是variant，每一行是sample
#The function eigen_windows() basically wants your data in a numeric matrix, with one row per variant and one column per sample (so that x[i,j] is the number of alleles that sample j has at site i). If your data are already in this form, then you can use it directly.
#2.Compute the distance matrix between windows - done with pc_dist() on the output of eigen_windows().
#3.Visualize - whatever you want; MDS is implemented in R with cmdscale().
#######首先是对VCF文件进行处理
#有两种模式，第一种是直接按照snps <- read_vcf("mydata.vcf")这样，把vcf文件直接读到内存里面
#还有一种，就是直接用bcftools转化成bcf并建立索引\
for chr in {1..42}
do 
	bcftools convert -O b chr${chr}_landEU.vcf.gz > chr${chr}_landEU.bcf &
done
for chr in {1..42}
do 
	bcftools index chr${chr}_landEU.bcf &
done
##########################***************************************************************************************************在服务器上面提交
#!/usr/bin/Rscript
# devtools::install_github("petrelharp/local_pca/lostruct")
library(tidyverse) #tidyverse包是对一些具有相同思想，且可以一同工作的R包的收集。
library(lostruct) ##这个是用来做local PCA的
library(colorspace)  ##这个是用来画图的
library(RColorBrewer) ##这个是用来画图的
library(ggmap)  ##这个是用来画地图的
library(zoo)
library(Matrix) ##这个是用来转换数据格式的
library(grid) ##这个是用来画地图的
library(scatterpie)
library(SNPRelate)	##
library(gridExtra)  ##这个是用来画地图的
library(ggExtra)  ##这个是用来画地图的
library(spdep)  ##这个包是用来做Rotation的
#开始读头文件，知道有哪些个体和那些位点 win100_region.r
for (chr in seq(1,42)){
	window_size <- 100 #这是指每100个snp算一个值
	k_kept <- 40  #这是指每40个window算一个mds
	max_distance_between_outliers <- 100
	samples <- read.table("/data2/xuebo/Projects/Speciation/group/landrace_group/landrace_EU.txt",header =F)
	colnames(samples) <- c("sample")
	bcf.file <- paste("/data2/xuebo/Projects/Speciation/lostruct/EU_land/vcffile/chr",chr, "_landEU.bcf",sep="")
	sites <- vcf_positions(bcf.file)
	#the function vcf_windower() will create the window extractor function 这个方法是按照窗口的大小提取数据
	win.fn.snp <- vcf_windower(bcf.file, size=window_size, type='snp', sites=sites) 
	win.fn.snp (3) #查看数据，这个是显示的是第三个window的情况
	system.time( snp.pca <- eigen_windows(win.fn.snp,k=2) )  ##system.time这个的作用是输出运行时间  #a matrix whose rows give the first two eigenvalues and eigenvectors for each window 每个window的特征值和特征向量
	system.time( pcdist <- pc_dist( snp.pca ) )   ##the pairwise distance matrix between those windows 这是生成的两两的值的矩阵
	na.inds <- is.na(pcdist[,1])   #统计方阵里面有多少个NA 
	  if (sum(na.inds) == length(na.inds)){
	    na.inds <- is.na(pcdist[,2]) 
	  }  #这个判断的含义是要是na.inds的个数和na.inds的长度相等，na.inds就是pcdist不是NA的第二列的值
	mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept )  #去掉方阵里面的NA并且每40个window算一个mds,mds是一个list，一共含有5列，points，eig，x，ac，GOF
	##现在是local PCA
	mds.coords <- mds$points #提取mds的其中一列，这一列是MDS coordinate
	colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords)) #加上名字
	win.regions <- region(win.fn.snp)() #这是得出的是SNP的region,即染色体和位置
	win.regions$n <- 1:nrow(win.regions) #每一列标上列数n，即set number,放在最后一列
	win.regions <- win.regions[!na.inds,]  #去掉NA
	win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions #加上一列mid,mutate()用于创建或修改列，现在的win.regions是这样的：chrom   start     end n       mid
	#Add the columns for all the MDS coordinates,下面这个小的循环是建立mds01-mds40,里面的值都是NA，其实就是加列
	  for (k in 1:k_kept){
	    str_pad(k, 2, pad = "0")   
	    name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
	    win.regions$tmp <- "NA"
	    win.regions <- win.regions %>% rename(!!name := tmp) 
	  }  
	#Add the MDS coordinates to each window.把刚才加上的40列填上数，填上的是mds.coords里面的值，每5个值取一个，放进mds01-mds40里面
	  for (i in 1:k_kept){
	    j = i + 5
	    win.regions[,j] <- mds.coords[,i]
	  }
	saveRDS(win.regions, file = paste("chr",chr,".",window_size,".windows.rds",sep=""))
}

#############################


































