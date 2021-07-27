`getHybrid_all` <- function(infile,times,out,head=F){
  dat = read.table(infile,head=head)
  dat[,1:ncol(dat)] = sapply(dat[,1:ncol(dat)],as.character)
  if(ncol(dat)>1){
    ng = unique(dat[,2])
  }else{
    ng = 1
  }
  
  ID1 = NULL
  ID2 = NULL
  for (i in (1:(nrow(dat)-1))){
    for (j in ((i+1):nrow(dat)) ){
      ID1 = append(ID1,dat[i,1])
      ID2 = append(ID2,dat[j,1])
    }
  }
  data = data.frame(ID1,ID2)
  write.table(data,paste(out,"_all",".txt",sep=""),col.names = F,row.names = F,quote=F,sep="\t")
}
`resample` <- function(infile,size,outfile){
  f = read.table(infile,head=F)
  # f= read.table(paste(outfile,"_all.txt",sep=""),header = F)
  mv = min(nrow(f),size)
  idx = sample(1:nrow(f),mv,replace = F)
  f1 = f[idx,]
  write.table(f1,paste(outfile,"_",mv,".txt",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
}
#Wild_einkorn
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Wild_einkorn.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Wild_einkorn"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
#Domesticated_einkorn
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Domesticated_einkorn.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Domesticated_einkorn"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
#Urartu
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Urartu.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Urartu"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
#speltoides
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Speltoides.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/speltoides"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
#Wild_emmer
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Wild_emmer.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/wild_emmer"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
#Domesticated_emmer
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Domesticated_emmer.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/dom_emmer"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Rivet_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Rivet_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Rivet_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Georgian_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Georgian_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Georgian_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Ispahanicum
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Ispahanicum.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Ispahanicum"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Polish_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Polish_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Polish_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Persian_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Persian_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Persian_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Khorasan_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Khorasan_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Khorasan_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Durum
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Durum.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Durum"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Spelt
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Spelt.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Spelt"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Macha
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Macha.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Macha"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Club_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Club_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Club_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Indian_dwarf_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Indian_dwarf_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/indian_dwarf"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Yunan
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Yunan_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/yunna"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Xinjiang_wheat
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Xinjiang_wheat.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Xinjiang_wheat"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Tibetan_semi_wild
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Tibetan_semi_wild.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Tibetan_semi_wild"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Vavilovii
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Vavilovii.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Vavilovii"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Landrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Landrace.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Landrace"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##EULandrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/landrace_group")
infile ="landrace_EU.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/EULandrace"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##WALandrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/landrace_group")
infile ="landrace_WA.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/WALandrace"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##SCALandrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/landrace_group")
infile ="landrace_SCA.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/SCALandrace"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##EALandrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/landrace_group")
infile ="landrace_EA.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/EALandrace"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Cultivar
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Cultivar.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Cultivar"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Strangulata
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Strangulata.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Strangulata"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Meyeri
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Meyeri.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Meyeri"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)
##Anathera
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/group/subspecies")
infile ="sub_Anathera.txt"
outfile = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/group2/hybrid_file2/Anathera"
getHybrid_all(infile,10,outfile)
resample(paste(outfile,"_all.txt",sep=""),100,outfile)


###冰川数据整理
library(ggplot2)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/Climate_change/")
data <- read.table("LR04stack.txt",header=T,stringsAsFactors = F,sep="\t")
ggplot(data, aes(x=data$Time..ka., y=data$Benthic.d18O..per.mil.,color='qsec')) +
  geom_line() +
  scale_x_log10()+scale_y_reverse() +
  geom_hline(aes(yintercept=4), colour="grey", linetype="dashed") +
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
  labs(x = "Time(Ka)",y = "O18 (proxy the temperature)")





#现在是对smc++的结果进行整理和统计
##install.packages("rjson")
dev.off()
library("rjson")
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV1")
json_data <- fromJSON(file="model.final.json")
x=c(0.065,
    0.08060878939574145,
    0.09996579888995383,
    0.1239711081461128,
    0.15374093765702604,
    0.19065955177075428,
    0.23644362546115683,
    0.27718075192442565,
    0.8461752431057717,
    1.4358649242509498,
    2.0478122329617787,
    2.6837635122081913,
    3.345679055025798,
    4.035769546957803,
    4.756540614051093,
    5.510847740950179,
    6.301964597926417,
    7.133668908223848,
    8.010351552870036,
    8.937156893179546,
    9.920165683644731,
    10.96663709577083,
    12.085334369569281,
    13.286971359425259,
    14.584838181489143,
    15.995699740514961,
    17.541123724389937,
    19.249510679333554,
    21.159324810108256,
    23.324493496032865,
    25.823997307953135,
    28.780282850853457,
    32.398483936571765,
    37.06315643441619,
    43.637643063035284,
    54.876802189498804,
    65.0)
y = c(2.0824050397778535,
      1.970063947492444,
      1.7654580396210888,
      1.5093472967621369,
      1.204778678355281,
      0.9052899833170523,
      0.7249759714446001,
      0.7434214484145799,
      -0.14722907912065877,
      -0.1210176904690863,
      0.34752210915433196,
      0.2926724485728082,
      0.06428106791216802,
      -0.17046149393317603,
      -0.37165616476738816,
      -0.5005080924708787,
      -0.5775714478372704,
      -0.6170495707637246,
      -0.6322667040062419,
      -0.6302567251057744,
      -0.6186835817973577,
      -0.6007587258659589,
      -0.5786514247918144,
      -0.5543009838274012,
      -0.5285608612890621,
      -0.5022097456132978,
      -0.4757662404210534,
      -0.4495788162389056,
      -0.42388289403374774,
      -0.3988409892113659,
      -0.37457228275734866,
      -0.3511816395788078,
      -0.32888813898286223,
      -0.30809147434299966,
      -0.28951221026478013,
      -0.2750880742833957,
      -0.2721106112231711)
N0 = 7692.307692307692
xp = 2*N0*x
yp = N0*exp(y)
plot(xp,yp, xlab="Generation", ylab="Ne", log="xy",type = "o",col = "#D4005F",pch=19,cex=1.6,lwd=2)
###
# sed  -i '' 5,64d  model.final.json 这是在mac上面的操作
library("rjson")
wildeinkorn <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/plot_para3_4/wildeinkorn/model.final.json")
N0_WEI = wildeinkorn$model$N0
xp_WEI = 2*N0*wildeinkorn$model$knots
yp_WEI = N0*exp(wildeinkorn$model$y)
plot(xp_WEI,yp_WEI, xlab="Generation", ylab="Ne", log="xy",type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,xlim = c(1000,1000000),ylim = c(1000,1000000))
domeinkorn <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/plot_para3_4/domeinkorn/model.final.json")
N0_DEI = domeinkorn$model$N0
xp_DEI = 2*N0*domeinkorn$model$knots
yp_DEI = N0*exp(domeinkorn$model$y)
points(xp_DEI,yp_DEI,type = "o",col = "blue",pch=19,cex=1,lwd=2)
urartu <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/plot_para3_4/urartu/model.final.json")
N0_urartu = urartu$model$N0
xp_urartu = 2*N0*urartu$model$knots
yp_urartu = N0*exp(urartu$model$y)
points(xp_urartu,yp_urartu,type = "o",col = "orange",pch=19,cex=1,lwd=2)
speltoides <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/plot_para3_4/speltoides/model.final.json")
N0_speltoides = speltoides$model$N0
xp_speltoides = 2*N0*speltoides$model$knots
yp_speltoides = N0*exp(speltoides$model$y)
points(xp_speltoides,yp_speltoides,type = "o",col = "green",pch=19,cex=1,lwd=2)
Strangulata <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/plot_para3_4/Strangulata/model.final.json")
N0_Strangulata = Strangulata$model$N0
xp_Strangulata = 2*N0*Strangulata$model$knots
yp_Strangulata = N0*exp(Strangulata$model$y)
points(xp_Strangulata,yp_Strangulata,type = "o",col = "pink",pch=19,cex=1,lwd=2)
cbind(xp_WEI,yp_WEI,xp_DEI,yp_DEI,xp_urartu,yp_urartu,xp_speltoides,yp_speltoides,xp_Strangulata,yp_Strangulata)


####现在是四个地区的landrace
library("rjson")
EU <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/EU_model.final.json")
N0_EU = EU$model$N0
xp_EU = 2*N0*EU$model$knots
yp_EU = N0*exp(EU$model$y)
plot(xp_EU,yp_EU, xlab="Generation", ylab="Ne", log="xy",type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,xlim = c(1000,200000),ylim = c(100,2000000))
boxplot(xp_EU)

WA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/WA_model.final.json")
N0_WA = WA$model$N0
xp_WA = 2*N0*WA$model$knots
yp_WA = N0*exp(WA$model$y)
points(xp_WA,yp_WA,type = "o",col = "orange",pch=19,cex=1,lwd=2)
SCA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/SCA_model.final.json")
N0_SCA = SCA$model$N0
xp_SCA = 2*N0*SCA$model$knots
yp_SCA = N0*exp(SCA$model$y)
points(xp_SCA,yp_SCA,type = "o",col = "blue",pch=19,cex=1,lwd=2)
EA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/EA_model.final.json")
N0_EA = EA$model$N0
xp_EA = 2*N0*EA$model$knots
yp_EA = N0*exp(EA$model$y)
points(xp_EA,yp_EA,type = "o",col = "pink",pch=19,cex=1,lwd=2)
allsize = cbind(xp_EU,yp_EU,xp_WA,yp_WA,xp_SCA,yp_SCA,xp_EA,yp_EA)

##开始画地图
library(tess3r)
library(RColorBrewer)
display.brewer.all()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace")
datasize= read.table("landrace_size.txt",header = T)
datasizetem1 = datasize$X1000_2000/(max(datasize$X1000_2000) + 10)
datasizetem11 = 1-datasizetem1
datasizetem2 = datasize$X2000_3000/(max(datasize$X2000_3000) + 10)
datasizetem22 = 1-datasizetem2
datasizetem3 = datasize$X3000_5000/(max(datasize$X3000_5000) + 10)
datasizetem33 = 1-datasizetem3
datasizetem4 = datasize$X5000_7000/(max(datasize$X5000_7000) + 10)
datasizetem44 = 1-datasizetem4
datasizetem5 = datasize$X7000_10000/(max(datasize$X7000_10000) + 10)
datasizetem55 = 1-datasizetem5

datasizeall = cbind(datasize,datasizetem1,datasizetem11,datasizetem2,datasizetem22,datasizetem3,datasizetem33,
                    datasizetem4,datasizetem44,datasizetem5,datasizetem55)
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/introgression/Dstatistic/individuals/1to1_indi/map/samples_figure1_noAM.csv",head =T,sep=",")
info = info1[!is.na(info1$latitude),]
index = match(datasizeall$Taxa,info$samples)
lati = info[index,7:8]
datasizemap = cbind(datasizeall,lati)
datasizemap = na.omit(datasizemap)
##
my.colors <- c("tomato", "#f5deda")
my.colors2 <- c("tomato", "#fa6c52")
my.palette <- CreatePalette(my.colors, 20)

my.colors3 = c(CreatePalette(my.colors, 20)[[1]])
my.palette2 <- list(my.colors3,CreatePalette(my.colors2, 40)[[1]][3:11])

#1
datamatrix_1= cbind(datasizemap$datasizetem1,datasizemap$datasizetem11)
datacoord_1 = cbind(datasizemap$longitude,datasizemap$latitude)
plot(as.qmatrix(datamatrix_1), datacoord_1, method = "map.max", interpol = FieldsKrigModel(50),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette2,main="1000-2000")
#2
datamatrix_2= cbind(datasizemap$datasizetem2,datasizemap$datasizetem22)
plot(as.qmatrix(datamatrix_2), datacoord_1, method = "map.max", interpol = FieldsKrigModel(50),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="2000-3000")
#3
datamatrix_3= cbind(datasizemap$datasizetem3,datasizemap$datasizetem33)
plot(as.qmatrix(datamatrix_3), datacoord_1, method = "map.max", interpol = FieldsKrigModel(50),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="3000-5000")
#4
datamatrix_4= cbind(datasizemap$datasizetem4,datasizemap$datasizetem44)
plot(as.qmatrix(datamatrix_4), datacoord_1, method = "map.max", interpol = FieldsKrigModel(50),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="5000-7000")
#5
datamatrix_5= cbind(datasizemap$datasizetem5,datasizemap$datasizetem55)
plot(as.qmatrix(datamatrix_5), datacoord_1, method = "map.max", interpol = FieldsKrigModel(50),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="7000-10000")


####现在是四个地区的landrace的D
library("rjson")
EU <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/EU_D_model.final.json")
N0_EU = EU$model$N0
xp_EU = 2*N0*EU$model$knots
yp_EU = N0*exp(EU$model$y)
plot(xp_EU,yp_EU, xlab="Generation", ylab="Ne", log="xy",type = "o",col = "#D4005F",pch=19,cex=1,lwd=2,xlim = c(1000,200000),ylim = c(50,2000000))
WA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/WA_D_model.final.json")
N0_WA = WA$model$N0
xp_WA = 2*N0*WA$model$knots
yp_WA = N0*exp(WA$model$y)
points(xp_WA,yp_WA,type = "o",col = "orange",pch=19,cex=1,lwd=2)
SCA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/SCA_D_model.final.json")
N0_SCA = SCA$model$N0
xp_SCA = 2*N0*SCA$model$knots
yp_SCA = N0*exp(SCA$model$y)
points(xp_SCA,yp_SCA,type = "o",col = "blue",pch=19,cex=1,lwd=2)
EA <- fromJSON(file="/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace/EA_D_model.final.json")
N0_EA = EA$model$N0
xp_EA = 2*N0*EA$model$knots
yp_EA = N0*exp(EA$model$y)
points(xp_EA,yp_EA,type = "o",col = "pink",pch=19,cex=1,lwd=2)
allsize = cbind(xp_EU,yp_EU,xp_WA,yp_WA,xp_SCA,yp_SCA,xp_EA,yp_EA)
##
library(tess3r)
library(RColorBrewer)
display.brewer.all()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/plotV2/test_landrace")
datasize= read.table("landrace_sizeD.txt",header = T)
datasizetem1 = datasize$X1000_2000/(max(datasize$X1000_2000) + 10)
datasizetem11 = 1-datasizetem1
datasizetem2 = datasize$X2000_3000/(max(datasize$X2000_3000) + 10)
datasizetem22 = 1-datasizetem2
datasizetem3 = datasize$X3000_5000/(max(datasize$X3000_5000) + 10)
datasizetem33 = 1-datasizetem3
datasizetem4 = datasize$X5000_7000/(max(datasize$X5000_7000) + 10)
datasizetem44 = 1-datasizetem4
datasizetem5 = datasize$X7000_10000/(max(datasize$X7000_10000) + 10)
datasizetem55 = 1-datasizetem5

datasizeall = cbind(datasize,datasizetem1,datasizetem11,datasizetem2,datasizetem22,datasizetem3,datasizetem33,
                    datasizetem4,datasizetem44,datasizetem5,datasizetem55)
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/introgression/Dstatistic/individuals/1to1_indi/map/samples_figure1_noAM.csv",head =T,sep=",")
info = info1[!is.na(info1$latitude),]
index = match(datasizeall$Taxa,info$samples)
lati = info[index,7:8]
datasizemap = cbind(datasizeall,lati)
datasizemap = na.omit(datasizemap)
##
my.colors <- c("f0c5c5", "#f5deda")
my.colors2 <- c("tomato", "#fa6c52")
my.palette <- CreatePalette(my.colors, 20)

my.colors3 = c(CreatePalette(my.colors, 20)[[1]])
my.palette2 <- list(my.colors3,CreatePalette(my.colors2, 40)[[1]][3:11])

#1
datamatrix_1= cbind(datasizemap$datasizetem1,datasizemap$datasizetem11)
datacoord_1 = cbind(datasizemap$longitude,datasizemap$latitude)
plot(as.qmatrix(datamatrix_1), datacoord_1, method = "map.max", interpol = FieldsKrigModel(100),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="1000-2000")
#2
datamatrix_2= cbind(datasizemap$datasizetem2,datasizemap$datasizetem22)
plot(as.qmatrix(datamatrix_2), datacoord_1, method = "map.max", interpol = FieldsKrigModel(100),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="2000-3000")
#3
datamatrix_3= cbind(datasizemap$datasizetem3,datasizemap$datasizetem33)
plot(as.qmatrix(datamatrix_3), datacoord_1, method = "map.max", interpol = FieldsKrigModel(100),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="3000-5000")
#4
datamatrix_4= cbind(datasizemap$datasizetem4,datasizemap$datasizetem44)
plot(as.qmatrix(datamatrix_4), datacoord_1, method = "map.max", interpol = FieldsKrigModel(100),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="5000-7000")
#5
datamatrix_5= cbind(datasizemap$datasizetem5,datasizemap$datasizetem55)
plot(as.qmatrix(datamatrix_5), datacoord_1, method = "map.max", interpol = FieldsKrigModel(100),
     xlab = "Longitude", ylab = "Latitude",resolution = c(500,500), cex =0.6,col.palette = my.palette,main="7000-10000")


##################现在做的是有效群体大小和Drift的关系
#https://stephens999.github.io/fiveMinuteStats/wright_fisher_model.html
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
#data.frame to be filled
wf_df <- data.frame()

# effective population sizes
sizes <- c(50, 100, 1000, 5000,10000,100000)

# starting allele frequencies
starting_p <- c(.01, .1, .5, .8)

# number of generations
n_gen <- 100

# number of replicates per simulation
n_reps <- 50

# run the simulations
for(N in sizes){
  for(p in starting_p){
    p0 <- p
    for(j in 1:n_gen){
      X <- rbinom(n_reps, 2*N, p)  #生成二项分布随机数的函数是：rbinom()
      p <- X / (2*N)
      rows <- data.frame(replicate = 1:n_reps, N = rep(N, n_reps), 
                         gen = rep(j, n_reps), p0 = rep(p0, n_reps), 
                         p = p)
      wf_df <- bind_rows(wf_df, rows)
    }
  }
}

# plot it up!
p <- ggplot(wf_df, aes(x = gen, y = p, group = replicate)) +
  geom_path(alpha = .5) + facet_grid(N ~ p0) + guides(colour=FALSE)
p

##########################现在做的是split-time的分析
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/split15")
Strang_Land = read.table("split_time_Strang_Landrace.txt")
boxplot(Strang_Land)
WA_EU = read.table("split_time_WA_EU.txt")
boxplot(WA_EU)
WA_EA = read.table("split_time_WA_EA.txt")
boxplot(WA_EA)
data <- data.frame(
  name=c( rep("Strang_Land",length(Strang_Land$V1)), rep("WA_EU",length(WA_EU$V1)), rep("WA_EA",length(WA_EA$V1))),
  value=c(Strang_Land$V1, WA_EU$V1, WA_EA$V1 )
)
sample_size = data %>% group_by(name) %>% summarize(num=n())
data1 = data %>%left_join(sample_size) %>%mutate(myaxis = paste0(name, "\n", "n=", num))
ggplot(data1,aes(x=myaxis, y=value, fill=name)) +
  geom_violin(width=1, color="grey", alpha=1) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")+theme_bw()


#######Polynomial curve fitting and confidence interval
# We create 2 vectors x and y. It is a polynomial function.
x <- runif(300, min=-30, max=30) 
y <- -1.2*x^3 + 1.1 * x^2 - x + 10 + rnorm(length(x),0,100*abs(x)) 

# Basic plot of x and y :
plot(x,y,col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3 , xlab="" , ylab="") 

# Can we find a polynome that fit this function ?
model <- lm(y ~ x + I(x^2) + I(x^3))

# I can get the features of this model :
#summary(model)
#model$coefficients
#summary(model)$adj.r.squared

#For each value of x, I can get the value of y estimated by the model, and the confidence interval around this value.
myPredict <- predict( model , interval="predict" )

#Finally, I can add it to the plot using the line and the polygon function with transparency.
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix , 1], col=2, lwd=2 )
polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

















