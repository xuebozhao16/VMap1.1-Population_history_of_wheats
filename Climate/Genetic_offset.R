##参考教程 https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-021-01526-9/MediaObjects/41559_2021_1526_MOESM1_ESM.pdf
###现在的代码是计算genetic offset
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(RColorBrewer)
library(rasterVis)
library(LEA)
library(adegenet)
library(maps)
display.brewer.all()
library(dismo)
library(gplots)
library(raster)
library(gdistance)
#####################首先是产生数据，改数据格式，搞个模板出来
library(raster)
library(geosphere)
library(tidyverse)
library(ggplot2)
library(igraph)
library(ggridges)
library(UpSetR)
library(here)
library(rasterVis)
setwd(here::here())
#source("/Users/xuebozhao/Downloads/spiritu-santi-Climate-Change-Genomics-ad0ab0f/datasets/code/climate_change_functions.R")
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/present")
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
writeRaster(temp1, filename="bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/wc2.1_10m_bio/wc2.1_10m_tmin_BCC-CSM2-MR_ssp126_2021-2040.tif")
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_30/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/wc2.1_10m_bio/wc2.1_10m_tmin_BCC-CSM2-MR_ssp126_2041-2060.tif")
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_50/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/wc2.1_10m_bio/wc2.1_10m_tmin_BCC-CSM2-MR_ssp126_2061-2080.tif")
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_70/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)
##
temp1 <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/wc2.1_10m_bio/wc2.1_10m_tmin_BCC-CSM2-MR_ssp126_2081-2100.tif")
writeRaster(temp1, filename="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_90/bio_1.asc", format = "ascii", datatype='INT4S', overwrite=TRUE)

#####################现在是开始读文件，跑Gradient Forest，看哪一个环境变量的权重比较大
# read the input, see above for details
gfData <- read.table("/Users/xuebozhao/Downloads/spiritu-santi-Climate-Change-Genomics-ad0ab0f/datasets/input/gradient_forest_input.csv",header =
                       T,sep="\t",row.names = "pop")
# first step separate data based on the category of SNPS. In the function described in section 9, we do it for only 1 set of SNP at the time.
candidate <- gfData[,grep("candSNPs",names(gfData))]
reference <- gfData [,grep("refSNPs",names(gfData))]
# create a table with the bioclimatic information of populations
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
# define the bioclimatic variables of interest
bioclimatic <- paste("bio_",1:19,sep = "")
# set the importance of the permutation distribution of each variable. Type help(gradientForest for more details)
maxLevel <- log2(0.368*nrow(candidate)/2)
# run the algorithm (gradient forest function) for each set of SNPs
if(FALSE){ # FALSE if there is no need to run the analyses #要是第一次做需要Run，但是要是gf_runs.R这个文件产生之后(10M),就不需要再来这么一通了
  gf_candidate <- gradientForest(cbind(present[,bioclimatic], candidate),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(candidate), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  gf_reference <- gradientForest(cbind(present[,bioclimatic], reference),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(reference), ntree=500, 
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  # combine the GF models into a list and save it
  gf_runs <- list(gf_reference=gf_reference,
                  gf_candidate=gf_candidate)
  if(!dir.exists("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output")){
    dir.create("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output")
  }
  save(gf_runs,file = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output/gf_runs.R")
}
load("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output/gf_runs.R")
gf_candidate <- gf_runs$gf_candidate
gf_reference <- gf_runs$gf_reference
# Once the gradient forest model has been constructed, the importance of each variables (variable contribution) to the model can be estimated.
# The following vector contains such information
bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)]
most_cand <- names(bio_cand[1])
barplot(bio_cand,las=2,cex.names=0.8,col=c(brewer.pal(4,'Set3'),rep("grey",15)),ylab="Weigthed importance (R-sqr)")

#####################现在是计算allele的替换情况
# We can extract the allele turnover as a function of a single predictor variable (in this case bio_9). 
# This can be done for the combined SNPs (Overall option) or individual SNPs (Species option). 
# Note: here we show allele turnover across the range of individual variables, but we provide options to explore allele turnover across all variables. 
temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand,
                            type=c("Overall"),standardize = T) # al candidate SNPs
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand,
                        type=c("Species"),standardize = T) #each individual candidate allele
temp_ref_overall <- cumimp(gf_reference,predictor = most_cand,
                           type=c("Overall"),standardize = T) #all neutral SNPs
temp_ref_SNP <- cumimp(gf_reference,predictor = most_cand,
                       type=c("Species"),standardize = T) #each individidual neutral SNPs

# the next code is run only to estimate the y axis limit so the final plot incorporates all data.
# This is just esthetics and to automatize the code
ylim <- NULL
for(j in 1:length(temp_cand_SNP)){ #test each SNP
  ylim <- c(ylim,max(temp_cand_SNP[[j]][[2]])) # get the maximum value for a SNP
}
#same for reference SNPs, since we will plot both types of SNPs, we add the new maximum values to those that were obtained for the candidate SNPs 
# (we do not create a new ylim object)
for(j in 1:length(temp_ref_SNP)){
  ylim <- c(ylim,max(temp_ref_SNP[[j]][[2]]))
}
ylim <- max(ylim)
# code to plot the overall and individual allele turnover functions across the candidate bio
par(mfrow=c(1,2))
par(mai=c(0.9,0.8,0.4,0))
plot(temp_cand_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),
     ylab="Cumulative importance",xlab= "bio_9")
for(j in 1:length(temp_cand_SNP)){
  lines(temp_cand_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_cand_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
par(mai=c(0.9,0.1,0.4,0.6),tcl=-0.2)
plot(temp_ref_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="",xlab= "bio_9", yaxt="n")
for(j in 1:length(temp_ref_SNP)){
  lines(temp_ref_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f = 0.9))
}
lines(temp_ref_overall,col=brewer.pal(5,'Set3')[4],lwd=4)

#####################现在是计算和环境关联的情况
# Populations can be classified based on environmental categories defined from the allele turnover functions across a given environmental gradient; 
#here we show the patterns across bio_9.
dev.off()
pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))])
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) # get the x (bio value) and y (predicted cumulative importance) values of each population
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand]))) # identify which populations grow above the mean
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) # identify which populations grow above the mean
# record the categories of populations (they are adapted to cold or warm conditions) for future analyses.
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm]) # create a list containing the name of populations that belong to the environmental clusters.
#This list will be used in later analyses to classify populations.
plot(temp_cand_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),
     ylab="Cumulative importance",xlab= paste("Most important variable(",most_cand,")",sep=""),main="Candidate SNPs")
#for each individual SNP add the line, orange and lightblue indicate adaptive and reference SNPs
for(j in 1:length(temp_cand_SNP)){
  lines(temp_cand_SNP[[j]],col=adjustcolor(brewer.pal(5,'Set3')[3],alpha.f =0.6))
}
lines(temp_cand_overall,col=brewer.pal(5,'Set3')[4],lwd=4)
# this time we add the points of the populations depending on whether they grow in the cold or warm adapted cluster. We plot with different colors
warm_col=rev(brewer.pal(5,'OrRd'))
cold_col=rev(brewer.pal(5,'Blues'))
id_c <- order(temp$bio[cold])
id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col))
id_w <- order(temp$bio[warm])
id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col))
points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg=rev(id_wcol),cex=1.5)
points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg=id_ccol,cex=1.5)

#####################现在是计算genetic offset了
mask <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/present/wc2.1_10m_bio_1.asc") %>% replace(.,.,0)
# We create matrices of present and future climate data that will be used to extrapolate the functions constructed with the Gradient Forest analysis across the geographic landscape. 
# Briefly, the function converts any given number of enviromental raster layers into a data frame.
turn_score <- data.frame(gfData[,c("X","Y",most_cand)],temp)
present_mat <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/present/")
future_mat_50 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_50/") 
# you can change the layers to have other periods. 
#Here it is 2050-RCP4.5. We used this so the change is not so drastic, and we can test the code.
future_mat_70 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_70/") 
#Here it is 2070-RCP8.5, that will be harsher.
#predict allelic data across the landscape
pred_paSNPs <- predict(gf_candidate,present_mat[grep("bio",names(present_mat))])
pred_paSNPs_future_50 <-
  predict(gf_candidate,future_mat_50[grep("bio",names(future_mat_50))])
pred_paSNPs_future_70 <-
  predict(gf_candidate,future_mat_70[grep("bio",names(future_mat_70))])
# finally we estimate the Euclidian distance between the two matrices, this is the genetic offset; the euclidian_distance function is also an accessory function
euclidian_50 <-euclidian_distance(proj_fut=pred_paSNPs_future_50,pred_pres=pred_paSNPs)
euclidian_70 <-euclidian_distance(proj_fut=pred_paSNPs_future_70,pred_pres=pred_paSNPs) 
# create a raster layer that contains the genetic offset of each pixel. We add this information in the mask raster created at the beginning
offset_ras_50 <- mask
offset_ras_50[present_mat$cell]<- euclidian_50 #the present_mat$cell contains the cell of each pixel in the distribution of the species
offset_ras_70 <- mask
offset_ras_70[present_mat$cell]<- euclidian_70
#obtain the genetic offset of the know populations
genetic_off_50 <- raster::extract(offset_ras_50,present[,1:2])
genetic_off_70 <- raster::extract(offset_ras_70,present[,1:2])
#create the tables that contain the coordinates of populations and the genetic offset
pop_vul_50 <- data.frame(present[,1:2], genetic_off_50)
pop_vul_70 <- data.frame(present[,1:2], genetic_off_70)
# Add the category of the genetic offset to a table for future analyses. 
pop_vul_50$temp <- NA
pop_vul_50$temp[warm]<-"warm"
pop_vul_50$temp[cold]<-"cold"
pop_vul_50$temp <- factor(pop_vul_50$temp,levels = c("cold","warm"))
pop_vul_70$temp <- NA
pop_vul_70$temp[warm]<-"warm"
pop_vul_70$temp[cold]<-"cold"
pop_vul_70$temp <- factor(pop_vul_70$temp,levels = c("cold","warm"))
# get max offset
max_val <- max(c(pop_vul_50$genetic_off_50,pop_vul_70$genetic_off_70))
# create a table that the genetic offset of each cell in the map
vulnerable_areas <- as.data.frame(offset_ras_50)
vulnerable_areas <- data.frame(cell=1:nrow(vulnerable_areas),vul=vulnerable_areas[,1])
#extract only the cells and values of areas where there is a predicted vulnerability
vulnerable_areas_50 <- vulnerable_areas[which(vulnerable_areas$vul>0),]
# same for 2070
vulnerable_areas <- as.data.frame(offset_ras_70)
vulnerable_areas <- data.frame(cell=1:nrow(vulnerable_areas),vul=vulnerable_areas[,1])
vulnerable_areas_70 <- vulnerable_areas[which(vulnerable_areas$vul>0),]
# create object that contains the different outputs of GF, these will be used for other analyses in the next sections
GO_objects <- list(pred_paSNPs=pred_paSNPs, #predicted genetic space in the present,and the two future (next two)
                   pred_paSNPs_future_50=pred_paSNPs_future_50,
                   pred_paSNPs_future_70=pred_paSNPs_future_70,
                   present_mat=present_mat, #environmental values of all the cells in the present
                   vulnerable_areas_50=vulnerable_areas_50, #cell that have a predicted genetic offset in the future (same for the next line)
                   vulnerable_areas_70=vulnerable_areas_70,
                   genetic_off_50=genetic_off_50, #genetic offset of populations in the present and the future
                   genetic_off_70=genetic_off_70,
                   pop_vul_50=pop_vul_50, #table containing containing coordinates and geneticoffset
                   pop_vul_70=pop_vul_70,
                   present=present, #table containing bioclim data of popualtions
                   temp_cand_overall=temp_cand_overall, # turnover functions
                   gf_candidate=gf_candidate, # gradient forest model
                   offset_ras_50=offset_ras_50, #raster with the offset
                   offset_ras_70=offset_ras_70,
                   gfData=gfData, #input data
                   bio_cand=most_cand, #bio climatic variable with the strongest contribution
                   categories=categories,
                   turn_score=turn_score) #categories of populations
# save the outputs from GF for futher analyses
save(GO_objects,file = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output/GO_objects.R")
#finally write the information to the conservation file
vul_gen_off <- cbind(pop_vul_50[c("X","Y","genetic_off_50")],pop_vul_70[c("genetic_off_70")])
# create dir if it does not exist
write.csv(vul_gen_off,file = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/output/vul_fg.csv",row.names = T) # create a data frame indicating all areas that have a genetic offset abvove the threshold
gen_off_stack <- stack(offset_ras_50,offset_ras_70)
names(gen_off_stack) <- paste(c("Year_2050","Year_2070"))
#plot genetic offset and the known populations according to the environmental cluster in which they grow
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),
                     xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                     main = "Genomic offset",par.settings=rasterVis::rasterTheme(rev(brewer.pal(5,'Blues'))))

##################现在开始展示genetic offset的Box
#plot the distribution of genetic offset according to the environment in which populations grow and estimate if they are significantly different
warm_col=brewer.pal(12,'Paired')[1]
cold_col=brewer.pal(12,'Paired')[2]
cold_col1=brewer.pal(12,'Paired')[3]
par(mfrow=c(1,2))
par(mai=c(0.6,0.8,0.6,0))
pop_vul_50$temp2 = as.factor(c(rep("SCA",6),rep("EA",3),rep("EU",6),rep("WA",5),rep("EA",3)))
plot(pop_vul_50$genetic_off_50 ~ pop_vul_50$temp2, ylab="Genetic offset", xlab=
       most_cand, col=brewer.pal(4,'Paired'),main="2050",ylim=c(0,max_val))
# now, let’s plot according to the North/South distribution. Again, here the order is increasing and positive, so the order is N->S
par(mai=c(0.6,0.1,0.6,0.6),tcl=-0.2)
#test the significance
fit <- aov(pop_vul_50$genetic_off ~pop_vul_50$temp)
p <- summary(fit)
#add legend
pop_vul_70$genetic_off_70_2 = pop_vul_70$genetic_off_7 + 0.00005
pop_vul_70$temp2 = as.factor(c(rep(c("SCA","WA","EU","EA"),5),rep("EU",3)))
plot(pop_vul_70$genetic_off_70_2 ~ pop_vul_70$temp2, ylab="Genetic offset", xlab=
       bio_cand,col=brewer.pal(4,'Paired'),main="2070",ylim=c(0,max_val),yaxt="n")
fit <- aov(pop_vul_70$genetic_off ~ pop_vul_70$temp)
p2 <- summary(fit)
























