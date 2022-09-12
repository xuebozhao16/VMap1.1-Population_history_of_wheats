##fsc的最佳模型的结果——D
######这是plot the likelihoods as boxplots
early_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/model/early/Dlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/model/ongoing/Dlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/model/diff/Dlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/model/recent/Dlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/model/no/Dlineage_PopDiv_no.lhoods")

par(mfrow=c(1,1))
library(RColorBrewer)
display.brewer.all()

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="D lineage")
        #main="Model comparison with Likelihood distributions in D lineage")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
#ggplot2上色,表示上色之后并不好看
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- data.frame(
  name=c( rep("Different",length(diff_geneflow)), rep("Recent",length(recent_geneflow)), rep("Early",length(early_geneflow)), 
          rep("Ongoing",length(ongoing_geneflow)), rep('No', length(no_geneflow))  ),
  value=c(diff_geneflow, recent_geneflow, early_geneflow,ongoing_geneflow, no_geneflow )
)
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Model comparison with Likelihood distributions") +
  xlab("Model comparison")

##fsc的最佳模型的结果——AB
######这是plot the likelihoods as boxplots
par(mfrow=c(2,3), oma=c(4.5,3, 0, 0), mar=c(2,1.6,2,1), cex=1)
##10
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/WE_DE_10")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_DE_10")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
##20
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/WE_FT_20")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_FT_20")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
##30
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/WE_LAN_30")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_LAN_30")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
##21
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/DE_FT_21")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage DE_FT_21")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
##31
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/DE_LAN_31")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage DE_LAN_31")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
##32
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/SFS_multi/FT_LAN_32")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,diff_geneflow,recent_geneflow,early_geneflow,ongoing_geneflow,
        no_geneflow, xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage FT_LAN_32")
axis(side=1,at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))


###现在看的的是WA和其他的关系
library(RColorBrewer)
display.brewer.all()
dev.off()
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/landrace225/fsc_landrace/lhoods")
par(mfrow=c(3,3), oma=c(4.5,3, 0, 0), mar=c(2,1.6,2,1), cex=1)
##10 ('WA', 'EU')
PREFIX="fit_10_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs EU landrace (EU)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##20 ('WA', 'CA')
PREFIX="fit_20_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs CA landrace (CA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##30 ('WA', 'SA')
PREFIX="fit_30_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs SA landrace (SA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##40 ('WA', 'EA_South_West')
PREFIX="fit_40_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs Southeast EA landrace (SE_EA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##50 ('WA', 'L1_EA_Northwest')
PREFIX="fit_50_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs Northwest EA landrace (NE_EA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##60 ('WA', 'L2_EA_Northcentral')
PREFIX="fit_60_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs Northcentral EA landrace (NC_EA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##70 ('WA', 'L3_EA_Northeast')
PREFIX="fit_70_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs Northeast EA landrace (NE_EA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##80 ('WA', 'L8_EA_Southeast')
PREFIX="fit_80_Dlineage_PopDiv_"
no_geneflow<-scan(paste(PREFIX,"no",".lhoods",sep = ""))
early_geneflow<-scan(paste(PREFIX,"early",".lhoods",sep = ""))
ongoing_geneflow<-scan(paste(PREFIX,"ongoing",".lhoods",sep = ""))
recent_geneflow<-scan(paste(PREFIX,"recent",".lhoods",sep = ""))
diff_geneflow<-scan(paste(PREFIX,"diff",".lhoods",sep = ""))

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="WA landrace (WA) vs Southeast EA landrace (SE_EA)")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))


##############################################################现在做的是DAF
##fsc的最佳模型的结果——D
######这是plot the likelihoods as boxplots
no_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_multi/Dlineage_PopDiv_no.lhoods")
early_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_multi/Dlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_multi/Dlineage_PopDiv_ongoing.lhoods")
recent_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_multi/Dlineage_PopDiv_recent.lhoods")
diff_geneflow<-scan("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_multi/Dlineage_PopDiv_diff.lhoods")
par(mfrow=c(1,1))
library(RColorBrewer)
display.brewer.all()

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
         xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="D lineage")
#main="Model comparison with Likelihood distributions in D lineage")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different")) #5*5


##fsc的最佳模型的结果——AB
######这是plot the likelihoods as boxplots
par(mfrow=c(2,3), oma=c(4.5,3, 0, 0), mar=c(2,1.6,2,1), cex=1)
##10
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/WE_DE_10")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
         xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_DE_10")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##20
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/WE_FT_20")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_FT_20")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##30
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/WE_LAN_30")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
        xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage WE_LAN_30")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##21
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/DE_FT_21")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
         xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage DE_FT_21")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##31
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/DE_LAN_31")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
         xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage DE_LAN_31")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))
##32
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_multi/FT_LAN_32")
early_geneflow<-scan("ABlineage_PopDiv_early.lhoods")
ongoing_geneflow<-scan("ABlineage_PopDiv_ongoing.lhoods")
diff_geneflow<-scan("ABlineage_PopDiv_diff.lhoods")
recent_geneflow<-scan("ABlineage_PopDiv_recent.lhoods")
no_geneflow<-scan("ABlineage_PopDiv_no.lhoods")

boxplot(range = 0,no_geneflow,early_geneflow,ongoing_geneflow,recent_geneflow,diff_geneflow,
         xlab="Model comparison", ylab="Likelihood",xaxt="n",col=brewer.pal(5,"Set2"),
        border=brewer.pal(5,"Set2"),main="AB lineage FT_LAN_32")
axis(side=1,at=1:5, labels=c("No","Early","Ongoing","Recent","Different"))



######################################################现在看一下fsc估计物种形成时间的情况
library(ggplot2)
library(dplyr)
library(viridis)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_early")
DAF_SFS_D = read.table("Time_run20.txt")
DAF_SFS_D$year = DAF_SFS_D$V1 + (14071-mean(DAF_SFS_D$V1)) 
data <- data.frame(
  name=c(rep("rep1",length(DAF_SFS_D$V1)), rep("rep2",length(DAF_SFS_D$year))),
  value=c(DAF_SFS_D$V1,DAF_SFS_D$year))
sample_size = data %>% group_by(name) %>% summarize(num=n())
data1 = data %>%left_join(sample_size) %>%mutate(myaxis = paste0(name, "\n", "n=", num))
data1$name <- factor(data1$name,levels=c("rep1","rep2"))
medians = aggregate(value ~  name, data1, median) ##现在是求medians
medians_25_75 = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
c = as.character(unique(data1$name))
for(i in c(1:length(c))){
  medians_25_75$value25[i] = quantile(data1[which(data1$name == c[i]),2],0.25)
  medians_25_75$value75[i] = quantile(data1[which(data1$name == c[i]),2],0.75)
}
means = aggregate(value ~  name, data1, FUN = "mean") ##现在是求95%CI
medians_95CI = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
cc = as.character(unique(data1$name))
for(i in c(1:length(cc))){
  aa = data1[which(data1$name == c[i]),2]
  medians_95CI$value25[i] = mean(aa) - qnorm(0.975)*sd(aa)/sqrt(20)
  medians_95CI$value75[i] = mean(aa) + qnorm(0.975)*sd(aa)/sqrt(20)
}

ggplot(data1,aes(x=name, y=value, fill=name)) +
  geom_violin(width=1.8, color="grey", alpha=0.9) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  #scale_y_continuous(limits = c(0, 10000))+
  scale_x_discrete(labels=c("rep1","rep2")) +
  ##medians
  stat_summary(fun=median, colour="orange", geom="point", shape=15, size=3, show.legend=FALSE) + 
  geom_text(data = medians, aes(label = floor(value), y = value-500),size=3.5,colour = "orange") + 
  ##medians_25
  geom_point(data = medians_25_75,aes(x=name, y=value25), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value25), y = value25 - 500),size=3.5,colour = "red") + 
  ##medians_75
  geom_point(data = medians_25_75,aes(x=name, y=value75), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value75), y = value75 + 500),size=3.5,colour = "red") + 
  ##95CI_left
  geom_point(data = medians_95CI,aes(x=name, y=value25), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value25), y = value25 - 200),size=3.5,colour = "blue") + 
  #95CI_right
  geom_point(data = medians_95CI,aes(x=name, y=value75), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value75), y = value75 + 200),size=3.5,colour = "blue") + 
  theme(
    legend.position="none",
    #plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major=element_line(colour=NA)
  ) +
  ggtitle("") + ylab("Split time (Year)") +
  xlab("") +
  theme_bw()+theme(panel.grid=element_blank())

######################################################现在看一下fsc估计渗入时间的情况——
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_early")
DAF_SFS_D = read.table("TimeD_introgression_run20.txt")
DAF_SFS_D$year = DAF_SFS_D$V1 + (9729-mean(DAF_SFS_D$V1)) 
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage/DAF_SFS_early_diff2")
DAF_SFS_AB_WE_LAN = read.table("TimeAB_WE_LAN_introgression_run20.txt")
DAF_SFS_AB_WE_LAN$year = DAF_SFS_AB_WE_LAN$V1 + (8919-mean(DAF_SFS_AB_WE_LAN$V1))
DAF_SFS_AB_DE_LAN = read.table("TimeAB_DE_LAN_introgression_run20.txt")
DAF_SFS_AB_DE_LAN$year = DAF_SFS_AB_DE_LAN$V1 + (7228-mean(DAF_SFS_AB_DE_LAN$V1))
DAF_SFS_AB_FT_LAN = read.table("TimeAB_FT_LAN_introgression_run20.txt")
DAF_SFS_AB_FT_LAN$year = DAF_SFS_AB_FT_LAN$V1 + (3203-mean(DAF_SFS_AB_FT_LAN$V1))
##
data <- data.frame(
  name=c(rep("DAF_SFS_D",length(DAF_SFS_D$year)), rep("DAF_SFS_AB_WE_LAN",length(DAF_SFS_AB_WE_LAN$year)),
         rep("DAF_SFS_AB_DE_LAN",length(DAF_SFS_AB_DE_LAN$year)),rep("DAF_SFS_AB_FT_LAN",length(DAF_SFS_AB_FT_LAN$year))),
  value=c(DAF_SFS_D$year,DAF_SFS_AB_WE_LAN$year,DAF_SFS_AB_DE_LAN$year,DAF_SFS_AB_FT_LAN$year))
sample_size = data %>% group_by(name) %>% summarize(num=n())
data1 = data %>%left_join(sample_size) %>%mutate(myaxis = paste0(name, "\n", "n=", num))
data1$name <- factor(data1$name,levels=c("DAF_SFS_D","DAF_SFS_AB_WE_LAN","DAF_SFS_AB_DE_LAN","DAF_SFS_AB_FT_LAN"))
medians = aggregate(value ~  name, data1, median) ##现在是求medians
medians_25_75 = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
c = as.character(unique(data1$name))
for(i in c(1:length(c))){
  medians_25_75$value25[i] = quantile(data1[which(data1$name == c[i]),2],0.25)
  medians_25_75$value75[i] = quantile(data1[which(data1$name == c[i]),2],0.75)
}
means = aggregate(value ~  name, data1, FUN = "mean") ##现在是求95%CI
medians_95CI = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
cc = as.character(unique(data1$name))
for(i in c(1:length(cc))){
  aa = data1[which(data1$name == c[i]),2]
  medians_95CI$value25[i] = mean(aa) - qnorm(0.975)*sd(aa)/sqrt(20)
  medians_95CI$value75[i] = mean(aa) + qnorm(0.975)*sd(aa)/sqrt(20)
}

ggplot(data1,aes(x=name, y=value, fill=name)) +
  geom_violin(width=1.4, color="grey", alpha=0.9) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  #scale_y_continuous(limits = c(0, 10000))+
  scale_x_discrete(labels=c("DAF_SFS_D","DAF_SFS_AB_WE_LAN","DAF_SFS_AB_DE_LAN","DAF_SFS_AB_FT_LAN")) +
  ##medians
  stat_summary(fun=median, colour="orange", geom="point", shape=15, size=3, show.legend=FALSE) + 
  geom_text(data = medians, aes(label = floor(value), y = value-500),size=3.5,colour = "orange") + 
  ##medians_25
  geom_point(data = medians_25_75,aes(x=name, y=value25), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value25), y = value25 - 500),size=3.5,colour = "red") + 
  ##medians_75
  geom_point(data = medians_25_75,aes(x=name, y=value75), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value75), y = value75 + 500),size=3.5,colour = "red") + 
  ##95CI_left
  geom_point(data = medians_95CI,aes(x=name, y=value25), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value25), y = value25 - 200),size=3.5,colour = "blue") + 
  #95CI_right
  geom_point(data = medians_95CI,aes(x=name, y=value75), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value75), y = value75 + 200),size=3.5,colour = "blue") + 
  theme(
    legend.position="none",
    #plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major=element_line(colour=NA)
  ) +
  ggtitle("") + ylab("Split time (Year)") +
  xlab("") +
  theme_bw()+theme(panel.grid=element_blank())


######################################################现在看一下fsc估计渗入时间的情况——EA landrace
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/Dlineage/DAF_SFS_early")
DAF_SFS_D = read.table("TimeD_introgression_run20.txt")
DAF_SFS_D$year = DAF_SFS_D$V1 + (9633-mean(DAF_SFS_D$V1)) 
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/fsc/ABlineage_EA")
DAF_SFS_AB_WE_LAN = read.table("TimeAB_WE_LAN_introgression_run20.txt")
DAF_SFS_AB_WE_LAN$year = DAF_SFS_AB_WE_LAN$V1 + (9015-mean(DAF_SFS_AB_WE_LAN$V1))
DAF_SFS_AB_DE_LAN = read.table("TimeAB_DE_LAN_introgression_run20.txt")
DAF_SFS_AB_DE_LAN$year = DAF_SFS_AB_DE_LAN$V1 + (7688-mean(DAF_SFS_AB_DE_LAN$V1))
DAF_SFS_AB_FT_LAN = read.table("TimeAB_FT_LAN_introgression_run20.txt")
DAF_SFS_AB_FT_LAN$year = DAF_SFS_AB_FT_LAN$V1 + (6174-mean(DAF_SFS_AB_FT_LAN$V1))
##
data <- data.frame(
  name=c(rep("DAF_SFS_D",length(DAF_SFS_D$year)), rep("DAF_SFS_AB_WE_LAN",length(DAF_SFS_AB_WE_LAN$year)),
         rep("DAF_SFS_AB_DE_LAN",length(DAF_SFS_AB_DE_LAN$year)),rep("DAF_SFS_AB_FT_LAN",length(DAF_SFS_AB_FT_LAN$year))),
  value=c(DAF_SFS_D$year,DAF_SFS_AB_WE_LAN$year,DAF_SFS_AB_DE_LAN$year,DAF_SFS_AB_FT_LAN$year))
sample_size = data %>% group_by(name) %>% summarize(num=n())
data1 = data %>%left_join(sample_size) %>%mutate(myaxis = paste0(name, "\n", "n=", num))
data1$name <- factor(data1$name,levels=c("DAF_SFS_D","DAF_SFS_AB_WE_LAN","DAF_SFS_AB_DE_LAN","DAF_SFS_AB_FT_LAN"))
medians = aggregate(value ~  name, data1, median) ##现在是求medians
medians_25_75 = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
c = as.character(unique(data1$name))
for(i in c(1:length(c))){
  medians_25_75$value25[i] = quantile(data1[which(data1$name == c[i]),2],0.25)
  medians_25_75$value75[i] = quantile(data1[which(data1$name == c[i]),2],0.75)
}
means = aggregate(value ~  name, data1, FUN = "mean") ##现在是求95%CI
medians_95CI = data.frame(name= as.character(unique(data1$name)), value25= c(1:2),value75= c(1:2))
cc = as.character(unique(data1$name))
for(i in c(1:length(cc))){
  aa = data1[which(data1$name == c[i]),2]
  medians_95CI$value25[i] = mean(aa) - qnorm(0.975)*sd(aa)/sqrt(20)
  medians_95CI$value75[i] = mean(aa) + qnorm(0.975)*sd(aa)/sqrt(20)
}

ggplot(data1,aes(x=name, y=value, fill=name)) +
  geom_violin(width=1.4, color="grey", alpha=0.9) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  #scale_y_continuous(limits = c(0, 10000))+
  scale_x_discrete(labels=c("DAF_SFS_D","DAF_SFS_AB_WE_LAN","DAF_SFS_AB_DE_LAN","DAF_SFS_AB_FT_LAN")) +
  ##medians
  stat_summary(fun=median, colour="orange", geom="point", shape=15, size=3, show.legend=FALSE) + 
  geom_text(data = medians, aes(label = floor(value), y = value-500),size=3.5,colour = "orange") + 
  ##medians_25
  geom_point(data = medians_25_75,aes(x=name, y=value25), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value25), y = value25 - 500),size=3.5,colour = "red") + 
  ##medians_75
  geom_point(data = medians_25_75,aes(x=name, y=value75), color = 'red',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_25_75, aes(label = floor(value75), y = value75 + 500),size=3.5,colour = "red") + 
  ##95CI_left
  geom_point(data = medians_95CI,aes(x=name, y=value25), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value25), y = value25 - 200),size=3.5,colour = "blue") + 
  #95CI_right
  geom_point(data = medians_95CI,aes(x=name, y=value75), color = 'blue',shape=16, size=2, show.legend=FALSE)+
  geom_text(data = medians_95CI, aes(label = floor(value75), y = value75 + 200),size=3.5,colour = "blue") + 
  theme(
    legend.position="none",
    #plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major=element_line(colour=NA)
  ) +
  ggtitle("") + ylab("Split time (Year)") +
  xlab("") +
  theme_bw()+theme(panel.grid=element_blank())



