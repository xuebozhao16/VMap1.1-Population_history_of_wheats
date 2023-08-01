#根据环境变量对样本进行聚类
library(cluster)
library(factoextra)
#聚类
#根据经纬度给样本聚类
data <- read.table("225env.txt", header=F,stringsAsFactors = F)
colname <- c("elevation","temp1","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9","temp10","temp11","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","Latitude","Logititude")
rownames(data) <- data[,1]
data2 <- data[!is.na(data$V6),-1]
data2 <- data2[which(data2$V23 > -40),]
colnames(data2) <- colname
#按列进行标准化并聚类
df = scale(data2,center = T,scale = T)
colnames(df) <- colname
#------------------------------kmeans聚类---------
#确定应该分几个cluster
data2<- data2[,c(1:22)]
mydata <- data2
mydata <- df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
#kmeans聚类，标准化后
data2<- data2[,c(1:22)]
km <- kmeans(data2, 13,iter.max = 10000) #用于画地图
fviz_cluster(km, data = df,
  #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
  ellipse.type = "euclid",
  star.plot = T, 
  repel = TRUE,
  ggtheme = theme_minimal()
)
data2$type <- km$cluster

------------------------计算类群环境变量---------
lat_mean <- tapply(data2[,21],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,22],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:219) {
  for(j in 1:8){
    if(data2[i,23] == j ){
      data2[i,24] <- as.numeric(lat_mean[j])
      data2[i,25] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"13_cluster.txt",sep="\t",row.names = T,quote=F)

----------------------------------------地图上展示聚类结果---------------------------------
#new <- data.frame(cluster1=km$centers[,21], cluster2=km$centers[,22],size = km$size)
data <- read.table("13_cluster.txt",header=T,stringsAsFactors = F)
library(maps)
library(ggplot2)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-50,80)
mp_13<-mp+geom_point(aes(x=data2$Logititude, y=data2$Latitude,color = as.factor(data2$type)))+
  scale_size(range=c(1,1))+ 
  theme_classic()
  
data <- read.table("type.txt",header=T,stringsAsFactors = F)
mp_13<- mp+geom_point(aes(x=data$cluster2, y=data$cluster1,size=data$Type))+
  #scale_size(range=c(1,1))+ 
  theme_classic()
mp_13