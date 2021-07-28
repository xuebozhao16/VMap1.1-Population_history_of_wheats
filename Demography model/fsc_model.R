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





library(ggplot2)
library(hrbrthemes)

# Create dummy data
data <- data.frame(
  cond = rep(c("condition_1", "condition_2"), each=10), 
  my_x = 1:100 + rnorm(100,sd=9), 
  my_y = 1:100 + rnorm(100,sd=16) 
)

# linear trend + confidence interval
ggplot(data, aes(x=my_x, y=my_y)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum() +
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))

























