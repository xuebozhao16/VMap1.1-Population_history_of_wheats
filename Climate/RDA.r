#计算RDA：环境变量和遗传变异
#工作目录：/data2/yafei/Project3/Vmap1.1/Out/VCF/VmapE6/Landrace/Select_taxa

#1:提取计算迁徙路径的样本的VCF文件
#vcftools --gzvcf D_Land.vcf.gz --keep Select_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > Select_taxa/D_Land_Select.vcf.gz

#2:提取基因区的VCF文件
#bedtools intersect -a Gene.gff3 -b A_Land_Select.vcf.gz -wb > A_Land_Select_gene.vcf
#cat VCF.header A_Land_Select_gene.vcf | bgzip -c > A_Land_Select_gene.vcf.gz

#3:合并vcf并随机选取3000个位点
#vcf-concat A_Land_Select_gene.vcf.gz B_Land_Select_gene.vcf.gz D_Land_Select_gene.vcf.gz | bgzip -c > All_gene.vcf.gz
#run_pipeline.pl -Xms10g -Xmx200g -vcf All_gene.vcf.gz -sortPositions -export All_gene.hmp.txt -exportType HapmapDiploid

#sed '1d' All.sort.hmp.txt | shuf -n 5000 > shuf_5000.hmp.txt
#head -n 1 All.sort.hmp.txt > shuf.header
#cat shuf.header shuf_5000.hmp.txt > shuf_5000.hmp.txt2
#mv shuf_5000.hmp.txt2 shuf_5000.hmp.txt
#run_pipeline.pl -SortGenotypeFilePlugin -inputFile shuf_5000.hmp.txt -outputFile shuf_5000.sort.hmp.txt -fileType Hapmap
#run_pipeline.pl -Xmx100g -fork1 -h shuf_5000.sort.hmp.txt -export -exportType VCF

#4:RDA分析
#输入文件:
#  env_table:数据框，行是样本名，列是环境变量
#genotype_table:数据框，行是样本名，列是位点位置

library(vegan)
setwd("/Users/guoyafei/Documents/个人项目/Project-2-Migration/migration/add_ZNdata/Environment/")
taxa <- read.table("select_taxa.txt",header=T,stringsAsFactors = F)
taxa_EA <- taxa[which(taxa$Region=="EA"),1]
taxa_WA <- taxa[which(taxa$Region=="WA"),1]
taxa_SCA <- taxa[which(taxa$Region=="SCA"),1]
taxa_AF <- taxa[which(taxa$Region=="AF"),1]
taxa_EU <- taxa[which(taxa$Region=="EU"),1]

phylum <- read.delim('All_noMiss_0.05_2000.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
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
rda_tb.scaling1 <- summary(rda_tb_all, scaling = 2)
rda_tb.scaling1

#若只关注局部环境数据，除了在原始表格中修改变量个数外，还可直接在 rda() 中指定
#rda_part <- rda(phylum~elevation+one+two+three+four+five+six+seven+eight+nine+ten+eleven+twelve+thirteen+fourteen+fifteen+sixteen+seventeen+eighteen+nineteen, data = env, scale = FALSE)

label <- read.table("select_taxa.txt", header=T, stringsAsFactors = F)
label$Region <- as.factor(label$Region)
label$cols = label$Region 
label$cols = gsub("EA","#97FFFF",label$cols)
label$cols = gsub("EU","#FFD700",label$cols)
label$cols = gsub("SCA","#FF6347",label$cols)
label$cols = gsub("WA","#8470FF",label$cols)
label$cols = gsub("AF","#D8BFD8",label$cols)
label$cols = gsub("AM","#838B8B",label$cols)
cols <- label$cols

plot(rda_tb_all, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
#points(rda_tb, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
points(rda_tb_all, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = scales::alpha(cols, 0.7), cex = 1)

legend("topright", c("AF", "AM", "EA","EU", "SCA", "WA"),pch=19,col = c("#D8BFD8","#838B8B","#97FFFF","#FFD700","#FF6347","#8470FF") )
text(rda_tb_all, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

phylum_EA <- phylum_hel[taxa_EA,]
phylum_WA <- phylum_hel[taxa_WA,]
phylum_SCA <- phylum_hel[taxa_SCA,]
phylum_EU <- phylum_hel[taxa_EU,]

temp_EA <- env_temp[taxa_EA,]
temp_WA <- env_temp[taxa_WA,]
temp_SCA <- env_temp[taxa_SCA,]
temp_EU <- env_temp[taxa_EU,]

prec_EA <- env_prec[taxa_EA,]
prec_WA <- env_prec[taxa_WA,]
prec_SCA <- env_prec[taxa_SCA,]
prec_EU <- env_prec[taxa_EU,]

env_EA <- env_all[taxa_EA,]
env_WA <- env_all[taxa_WA,]
env_SCA <- env_all[taxa_SCA,]
env_EU <- env_all[taxa_EU,]

EU_temp_rda <- rda(phylum_EU~., temp_EU, scale = FALSE)
EU_temp.scaling1 <- summary(EU_temp_rda, scaling = 2)
EU_temp.scaling1
RsquareAdj(EU_temp_rda)

EU_prec_rda <- rda(phylum_EU~., prec_EU, scale = FALSE)
EU_prec.scaling1 <- summary(EU_prec_rda, scaling = 2)
EU_prec.scaling1
RsquareAdj(EU_prec_rda)

EU_all_rda <- rda(phylum_EU~., env_EU, scale = FALSE)
EU_all.scaling1 <- summary(EU_all_rda, scaling = 2)
EU_all.scaling1
RsquareAdj(EU_all_rda)

taxa_EA_25 <- c("TW030","TW032","ZN160","TW033","TW154","TW166","TW138","ZN177","ZN001","TW153","TW158","TW167","ZN083","TW160","TW155","TW151","TW159","TW162","TW145","TW144","ZN167","ZN097","ZN166","TW134","TW147")
taxa_SCA_25 <- c("ZN112","ZN179","TW029","TW051","TW025","TW028","TW094","TW027","ZN111","ZN119","TW057","XI_33","TW074","TW073","TW056","TW001","TW055","TW102","ZN113","ZN118","TW113","ZN121","XI_36","ZN115","TW087")
taxa_EU_25 <- c("TW085","XI_10","TW070","XI_7","TW065","XI_8","TW089","XI_5","XI_6","XI_4","TW052","TW071","TW059","TW108","TW062","TW072","TW026","XI_1","TW058","TW096","TW107","TW064","TW091","TW097","TW099")
taxa_WA_25 <- c("ZN176","TW078","ZN175","XI_31","XI_32","TW076","TW100","TW077","TW101","TW079","TW080","TW081","TW082","TW104","TW075","TW103","TW105","ZN110","ZN174","TW002","TW003","TW106","XI_16","XI_17","TW069")

phylum_EA <- phylum_hel[taxa_EA_25,]
phylum_WA <- phylum_hel[taxa_WA_25,]
phylum_SCA <- phylum_hel[taxa_SCA_25,]
phylum_EU <- phylum_hel[taxa_EU_25,]

temp_EA <- env_temp[taxa_EA_25,]
temp_WA <- env_temp[taxa_WA_25,]
temp_SCA <- env_temp[taxa_SCA_25,]
temp_EU <- env_temp[taxa_EU_25,]

prec_EA <- env_prec[taxa_EA_25,]
prec_WA <- env_prec[taxa_WA_25,]
prec_SCA <- env_prec[taxa_SCA_25,]
prec_EU <- env_prec[taxa_EU_25,]

env_EA <- env_all[taxa_EA_25,]
env_WA <- env_all[taxa_WA_25,]
env_SCA <- env_all[taxa_SCA_25,]
env_EU <- env_all[taxa_EU_25,]

#EU
EU_temp_rda <- rda(phylum_EU~., temp_EU, scale = FALSE)
RsquareAdj(EU_temp_rda)

EU_prec_rda <- rda(phylum_EU~., prec_EU, scale = FALSE)
RsquareAdj(EU_prec_rda)

EU_all_rda <- rda(phylum_EU~., env_EU, scale = FALSE)
RsquareAdj(EU_all_rda)

#EA
EA_temp_rda <- rda(phylum_EA~., temp_EA, scale = FALSE)
RsquareAdj(EA_temp_rda)

EA_prec_rda <- rda(phylum_EA~., prec_EA, scale = FALSE)
RsquareAdj(EA_prec_rda)

EA_all_rda <- rda(phylum_EA~., env_EA, scale = FALSE)
RsquareAdj(EA_all_rda)

#SCA
SCA_temp_rda <- rda(phylum_SCA~., temp_SCA, scale = FALSE)
RsquareAdj(SCA_temp_rda)

SCA_prec_rda <- rda(phylum_SCA~., prec_SCA, scale = FALSE)
RsquareAdj(SCA_prec_rda)

SCA_all_rda <- rda(phylum_SCA~., env_SCA, scale = FALSE)
RsquareAdj(SCA_all_rda)

#WA
WA_temp_rda <- rda(phylum_WA~., temp_WA, scale = FALSE)
RsquareAdj(WA_temp_rda)

WA_prec_rda <- rda(phylum_WA~., prec_WA, scale = FALSE)
RsquareAdj(WA_prec_rda)

WA_all_rda <- rda(phylum_WA~., env_WA, scale = FALSE)
RsquareAdj(WA_all_rda)

taxa_north <- c("TW150","TW001","ZN167","TW135","ZN097","TW149","TW148","TW161","ZN166","TW142","TW102","TW139","ZN113","TW113","TW173","TW086","TW141","TW147","TW168","TW132","XI_36","ZN155","TW165","TW087","TW137")
taxa_south <- c("TW029","TW094","TW028","ZN179","TW025","ZN177","TW051","ZN001","TW030","TW033","TW034","TW166","TW032","ZN160","TW167","ZN081","TW160","TW163","TW138","TW158","TW156","TW154","TW153","TW157")
taxa_south_west <- c("TW029","TW094","TW028","XI_15","TW025","ZN178","TW051","ZN069","TW031","TW033","TW034","TW172","TW032","ZN164","TW167")

phylum_north <- phylum_hel[taxa_north,]
phylum_south <- phylum_hel[taxa_south,]
phylum_south_west <- phylum_hel[taxa_south_west,]

temp_north <- env_temp[taxa_north,]
temp_south <- env_temp[taxa_south,]
temp_south_west <- env_temp[taxa_south_west,]

prec_north <- env_prec[taxa_north,]
prec_south <- env_prec[taxa_south,]
prec_south_west <- env_prec[taxa_south_west,]

env_north <- env_all[taxa_north,]
env_south <- env_all[taxa_south,]
env_south_west <- env_all[taxa_south_west,]

#north
north_temp_rda <- rda(phylum_north~., temp_north, scale = FALSE)
RsquareAdj(north_temp_rda)

north_prec_rda <- rda(phylum_north~., prec_north, scale = FALSE)
RsquareAdj(north_prec_rda)

north_all_rda <- rda(phylum_north~., env_north, scale = FALSE)
RsquareAdj(north_all_rda)

#South-west
South_west_temp_rda <- rda(phylum_south_west~., temp_south_west, scale = FALSE)
RsquareAdj(South_west_temp_rda)

South_west_prec_rda <- rda(phylum_south_west~., prec_south_west, scale = FALSE)
RsquareAdj(South_west_prec_rda)

#画各区域RDA分解图
plot(EU_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(EU_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(EU_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(WA_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(WA_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(WA_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(SCA_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(SCA_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(SCA_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(north_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(north_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(north_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(South_west_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(South_west_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(South_west_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)




#画barplot
data <- read.table("plot_data.txt", header=T, row.names = 1,stringsAsFactors = F)

barplot(data, beside = TRUE,
        col = c("lightblue", "mistyrose"),
        legend = rownames(data), ylim = c(0, 0.3))

