###这是下数据的网站 https://www.worldclim.org/data/worldclim21.html
library(raster)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/wc2.1_2.5m_bio")
temp1 <- raster("wc2.1_2.5m_bio_1.tif")
temp2 <- raster("wc2.1_2.5m_bio_2.tif")
temp3 <- raster("wc2.1_2.5m_bio_3.tif")
temp4 <- raster("wc2.1_2.5m_bio_4.tif")
temp5 <- raster("wc2.1_2.5m_bio_5.tif")
temp6 <- raster("wc2.1_2.5m_bio_6.tif")
temp7 <- raster("wc2.1_2.5m_bio_7.tif")
temp8 <- raster("wc2.1_2.5m_bio_8.tif")
temp9 <- raster("wc2.1_2.5m_bio_9.tif")
temp10 <- raster("wc2.1_2.5m_bio_10.tif")
temp11 <- raster("wc2.1_2.5m_bio_11.tif")
temp12 <- raster("wc2.1_2.5m_bio_12.tif")
temp13 <- raster("wc2.1_2.5m_bio_13.tif")
temp14 <- raster("wc2.1_2.5m_bio_14.tif")
temp15 <- raster("wc2.1_2.5m_bio_15.tif")
temp16 <- raster("wc2.1_2.5m_bio_16.tif")
temp17 <- raster("wc2.1_2.5m_bio_17.tif")
temp18 <- raster("wc2.1_2.5m_bio_18.tif")
temp19 <- raster("wc2.1_2.5m_bio_19.tif")
location <- read.csv("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/accessions_414.csv",header=T,sep=",",stringsAsFactors = F)
location2 <- location[1:414,]
site <- location2$Origin_country
lon <- location2$Logititude
lat <- location2$Latitude
samples <- data.frame(lon, lat)
temp.data <- samples 
temp.data$BIO1 <- extract(temp1, samples)
temp.data$BIO2 <- extract(temp2, samples)
temp.data$BIO3 <- extract(temp3, samples)
temp.data$BIO4 <- extract(temp4, samples)
temp.data$BIO5 <- extract(temp5, samples)
temp.data$BIO6 <- extract(temp6, samples)
temp.data$BIO7 <- extract(temp7, samples)
temp.data$BIO8 <- extract(temp8, samples)
temp.data$BIO9 <- extract(temp9, samples)
temp.data$BIO10 <- extract(temp10, samples)
temp.data$BIO11 <- extract(temp11, samples)
temp.data$BIO12 <- extract(temp12, samples)
temp.data$BIO13 <- extract(temp13, samples)
temp.data$BIO14 <- extract(temp14, samples)
temp.data$BIO15 <- extract(temp15, samples)
temp.data$BIO16 <- extract(temp16, samples)
temp.data$BIO17 <- extract(temp17, samples)
temp.data$BIO18 <- extract(temp18, samples)
temp.data$BIO19 <- extract(temp19, samples)
head(temp.data)
write.csv(temp.data, "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/414_bioclim.csv")

###现在是为了得到全球的数据点
a = seq(-180,180,2)
b = seq(-90,90,2)
ra = rep(a,91)
c = rep(1,length(ra))
for( i in 1:91){
  c[((i-1)*181+1):(i*181)]= b[i]
}
geo = data.frame(lat=ra,lon=c)
##
a = seq(-180,180,4)
b = seq(-90,90,2)
ra = rep(a,91)
c = rep(1,length(ra))
for( i in 1:91){
  c[((i-1)*91+1):(i*91)]= b[i]
}
geo = data.frame(lat=ra,lon=c)
###
a = seq(-180,180,8)
b = seq(-90,90,4)
ra = rep(a,46)
c = rep(1,length(ra))
for( i in 1:46){
  c[((i-1)*46+1):(i*46)]= b[i]
}
geo = data.frame(lat=ra,lon=c)
library(raster)
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/wc2.1_2.5m_bio")
temp1 <- raster("wc2.1_2.5m_bio_1.tif")
temp2 <- raster("wc2.1_2.5m_bio_2.tif")
temp3 <- raster("wc2.1_2.5m_bio_3.tif")
temp4 <- raster("wc2.1_2.5m_bio_4.tif")
temp5 <- raster("wc2.1_2.5m_bio_5.tif")
temp6 <- raster("wc2.1_2.5m_bio_6.tif")
temp7 <- raster("wc2.1_2.5m_bio_7.tif")
temp8 <- raster("wc2.1_2.5m_bio_8.tif")
temp9 <- raster("wc2.1_2.5m_bio_9.tif")
temp10 <- raster("wc2.1_2.5m_bio_10.tif")
temp11 <- raster("wc2.1_2.5m_bio_11.tif")
temp12 <- raster("wc2.1_2.5m_bio_12.tif")
temp13 <- raster("wc2.1_2.5m_bio_13.tif")
temp14 <- raster("wc2.1_2.5m_bio_14.tif")
temp15 <- raster("wc2.1_2.5m_bio_15.tif")
temp16 <- raster("wc2.1_2.5m_bio_16.tif")
temp17 <- raster("wc2.1_2.5m_bio_17.tif")
temp18 <- raster("wc2.1_2.5m_bio_18.tif")
temp19 <- raster("wc2.1_2.5m_bio_19.tif")
lon <- geo$lon
lat <- geo$lat
samples <- data.frame(lon, lat)
temp.data <- samples 
temp.data$BIO1 <- extract(temp1, samples)
temp.data$BIO2 <- extract(temp2, samples)
temp.data$BIO3 <- extract(temp3, samples)
temp.data$BIO4 <- extract(temp4, samples)
temp.data$BIO5 <- extract(temp5, samples)
temp.data$BIO6 <- extract(temp6, samples)
temp.data$BIO7 <- extract(temp7, samples)
temp.data$BIO8 <- extract(temp8, samples)
temp.data$BIO9 <- extract(temp9, samples)
temp.data$BIO10 <- extract(temp10, samples)
temp.data$BIO11 <- extract(temp11, samples)
temp.data$BIO12 <- extract(temp12, samples)
temp.data$BIO13 <- extract(temp13, samples)
temp.data$BIO14 <- extract(temp14, samples)
temp.data$BIO15 <- extract(temp15, samples)
temp.data$BIO16 <- extract(temp16, samples)
temp.data$BIO17 <- extract(temp17, samples)
temp.data$BIO18 <- extract(temp18, samples)
temp.data$BIO19 <- extract(temp19, samples)
head(temp.data)
dataall = na.omit(temp.data)
write.csv(dataall, "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/180X90/worldmap_bioclim3.csv")


####现在做的是Environmental_PCA
dev.off()
library(RColorBrewer)
display.brewer.all()
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/PCA/treeRangeColor2.txt",head =F)
brewer.pal(8,'Set1')
brewer.pal(11,'Spectral')
cols = c("#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3", "#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF",
         "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4" ,
         "#66C2A5" ,"#3288BD", "#5E4FA2")
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim")
dat = read.delim("accessions_414_bioclim.csv",sep=",",head=T)
data = dat[,c(1,12:30)]
dt = na.omit(data)
pca = prcomp(dt[,-1])
## PC
PCs = pca$x
### variance explained
s = summary(pca)
vrianceExplained = s$importance[2,]
## Plot
dim(dt)
dim(PCs)
index = match(pca_color$V1,dt$ID)
plot(PCs[index,1:2],col= "#BEBEBE",pch=19,bg="transparent",cex=1.5,
     xlab= paste("PC1 (", round(vrianceExplained[1]*100,2),"% )",sep=""),
     ylab= paste("PC2 (", round(vrianceExplained[2]*100,2),"% )",sep=""))
land = read.delim("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/group/group_fromYao/group_landraces.txt",head=F)
index2 = match(pca_color$V1,land$V1)
points(PCs[index2,1:2],col= "red",pch=19,cex=1.5)
cul = read.delim("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/group/group_fromYao/group_cultivar.txt",head=F)
index3 = match(pca_color$V1,cul$V1)
points(PCs[index3,1:2],col= "#0000FF",pch=19,cex=1.5)

##去掉landrace和cultivar再算
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim")
dat = read.delim("accessions_nobreadwheat_bioclim.csv",sep=",",head=T)
data = dat[,c(1,12:30)]
dt = na.omit(data)
pca = prcomp(dt[,-1])
## PC
PCs = pca$x
head(PCs)
dim(dat)
### variance explained
s = summary(pca)
vrianceExplained = s$importance[2,]
## Plot
index = match(pca_color$V1,dt$ID)
plot(PCs[index,1:2],col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5,
     xlab= paste("PC1 (", round(vrianceExplained[1]*100,2),"% )",sep=""),
     ylab= paste("PC2 (", round(vrianceExplained[2]*100,2),"% )",sep=""))



#############################################################开始画地图
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim")
dat = read.csv("accessions_414_bioclim.csv",head=T)
dat$Logititude = as.numeric(as.character(dat$Logititude))
data = dat[,c(1,7,8,12:30)]
dt = na.omit(data)
pca = prcomp(dt[,-(1:3)])
## PC
PCs = pca$x
### variance explained
s = summary(pca)
vrianceExplained = s$importance[2,]
## Plot
dim(dt)
dim(PCs)
index = match(pca_color$V1,dt$ID)
plot(PCs[index,1:2],col= "#BEBEBE",pch=19,bg="transparent",cex=1.5,
     xlab= paste("PC1 (", round(vrianceExplained[1]*100,2),"% )",sep=""),
     ylab= paste("PC2 (", round(vrianceExplained[2]*100,2),"% )",sep=""))
apply(pca$x, 1,sum)
pcasum = apply(abs(pca$x[,1:2]), 1,sum)
pcaper = abs(pca$x[,1:2])/pcasum
##
pca11 = pca$x[,1] + (abs(min(pca$x[,1])))
pca22  = pca11/max(pca11)
pcamatrix = cbind(pca22,1-pca22)
###
datacoord = cbind(dt$Logititude,dt$Latitude)
library(tess3r)
## Spatial interpolation of ancestry coefficient
my.colors <- c("orange","tomato")
my.palette <- CreatePalette(my.colors, 30)
barplot(t(as.matrix(pcamatrix)), border = NA, space = 0,xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix",col = my.colors) 
plot(as.qmatrix(pcamatrix), datacoord, method = "map.max",interpol = FieldsKrigModel(30),
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(500,500), cex = .4,
     col.palette = my.palette)

############
library(fields)
library(RColorBrewer)
library(mapplots)
library (LEA)
library(raster)
library(dplyr)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

###############
### Load data ###
###############

allData =read.table("path_to_your_file\your_file.txt", sep = "\t", header = TRUE)
# Required fields: latitude, longitude, PCoA1 and PCoA2 #
asc.raster='/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/paux/world.asc'
#pcoaData <- filter(allData, PCoA1 !="NA")
# Filtering out missing data #

long.all = dt$Logititude
lat.all = dt$Latitude
pcoa1.all =  PCs[,1]

asc.raster='world.asc'
# See attached archive for this file #

####################
### Perform analysis ###
####################

grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)
coordinates = cbind.data.frame(long=long.all,lat=lat.all)
cell_value_min = NULL
cell_value_max = NULL
colpalette = colorRampPalette(c("blue","cyan","orange","red","darkred"))(100)
mapadd=T
pts.size = .2
pts.shape = 19
theta=10
cluster = data.frame(pcoa1.all)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
###############
### Draw map ###
###############
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(pcoa1.all),max(pcoa1.all)))
# You can change de xlim and ylim to focus on a given region #
points(coordinates,pch=pts.shape,cex=0.1)
# Plot accessions coordinates on the map #



#####使用全球的数据画一次
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/180X90")
dat = read.csv("worldmap_bioclim3.csv",head=T)
dat$lon = as.numeric(as.character(dat$lon))
dat$lat = as.numeric(as.character(dat$lat))
dt = na.omit(dat)
pca = prcomp(dt[,-(1:3)])
## PC
PCs = pca$x
### variance explained
s = summary(pca)
vrianceExplained = s$importance[2,]
## Plot
dim(dt)
dim(PCs)
index = match(pca_color$V1,dt$ID)
plot(PCs[,1:2],col= "#BEBEBE",pch=19,bg="transparent",cex=1.5,
     xlab= paste("PC1 (", round(vrianceExplained[1]*100,2),"% )",sep=""),
     ylab= paste("PC2 (", round(vrianceExplained[2]*100,2),"% )",sep=""))
##
c = rep(1,18)
for(i in 4:22){
  c[i-3] = cor(PCs[,1],dt[,i])
}
max(c)


##
library(fields)
library(RColorBrewer)
library(mapplots)
library (LEA)
library(raster)
library(dplyr)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
asc.raster='/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim/paux/world.asc'
long.all = dt$lon
lat.all = dt$lat
pcoa1.all =  PCs[,1]
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)
coordinates = cbind.data.frame(long=long.all,lat=lat.all)
cell_value_min = NULL
cell_value_max = NULL
colpalette = colorRampPalette(c("blue","cyan","orange","red","darkred"))(100)
mapadd=T
pts.size = .2
pts.shape = 19
theta=10
cluster = data.frame(pcoa1.all)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(pcoa1.all),max(pcoa1.all)))
#414
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/worldclim")
dat414 = read.csv("accessions_414_bioclim.csv",head=T)
dat414$Logititude = as.numeric(as.character(dat414$Logititude))
data414 = dat414[,c(1,7,8,12:30)]
dt414 = na.omit(data414)
coord414 = cbind.data.frame(long=dt414$Logititude,lat=dt414$Latitude)
points(coord414,pch=pts.shape,cex=0.2) #12*12

##纯属好奇看一下气温和降水
cluster = data.frame(dt$BIO1)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(dt$BIO1),max(dt$BIO1)))
cluster = data.frame(dt$BIO5)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(dt$BIO5),max(dt$BIO5)))
cluster = data.frame(dt$BIO12)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(dt$BIO12),max(dt$BIO12)))
cluster = data.frame(dt$BIO16)
fit = Krig(coordinates,cluster,m = 1,theta = theta)
look<- predict(fit,grid)
out<- as.surface( grid, look)
if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
ncolors=length(colpalette)
plot(raster(out),col=colpalette,axes=F,xlim=c(-180,180),ylim=c(-90,90), zlim=c(min(dt$BIO16),max(dt$BIO16)))









