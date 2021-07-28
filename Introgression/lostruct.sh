###################Local PCA
######使用测试数据，找的是EUlanarace,2A,为了数据量少一些，使用的是带着外类群的/data2/xuebo/Projects/Speciation/tree/withBarley_segregate
bcftools view -Oz -S /data2/xuebo/Projects/Speciation/group/landrace_group/landrace_EU.txt /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr7.withBarley.vcf.gz --threads 48 -o chr7_landEU.withBarley.vcf.gz
for chr in {1..42}
do
	bcftools view -Oz -S /data2/xuebo/Projects/Speciation/group/landrace_group/landrace_EU.txt /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz --threads 48 -o chr${chr}_landEU.vcf.gz &
done
###
bcftools view  -Oz chr7_landEU.withBarley.vcf.gz -o chr7_landEU.withBarley.bcf.gz
# devtools::install_github("petrelharp/local_pca/lostruct")
library(tidyverse)
library(lostruct)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(zoo)
library(Matrix)
library(grid)
library(scatterpie)
library(SNPRelate)
library(gridExtra)
library(ggExtra)
samples <- read.table("/data2/xuebo/Projects/Speciation/lostruct/landrace_EU.txt",header =F)
colnames(samples) <- c("sample")
sites <- vcf_positions("/data2/xuebo/Projects/Speciation/lostruct/try/chr7_landEU.withBarley.bcf.gz")
win.fn.snp <- vcf_windower("/data2/xuebo/Projects/Speciation/lostruct/try/chr7_landEU.withBarley.bcf.gz", size=100, type="snp", sites=sites) 
#snps <- read_vcf("chr7_landEU.withBarley.vcf") #这是在跑bcftools
#eigenstuff <- eigen_windows(snps, win=100, k=2)  #这一步是多线程
system.time( snp.pca <- eigen_windows(win.fn.snp,k=2) )  ##mc.cores=这是多少个线程
system.time( pcdist <- pc_dist(snp.pca) )
#windist <- pc_dist(eigenstuff, npc=2 )
pcdist_na <- which(is.na(pcdist), TRUE)
na.inds <- is.na(pcdist[,1]) 
  if (sum(na.inds) == length(na.inds)){
    na.inds <- is.na(pcdist[,2]) 
  }
mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=40 )
##现在是整体的PCA
pdf("chr7_landEU.MDSplots1.pdf", height=6,width=6)
plot( mds$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(windist)) )
dev.off()
##现在是local PCA
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions$n <- 1:nrow(win.regions)
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
#Add the columns for all the MDS coordinates
  for (k in 1:40){
    str_pad(k, 2, pad = "0")
    
    name = paste("mds",str_pad(k, 2, pad = "0"),sep="");s
    win.regions$tmp <- "NA"
    win.regions <- win.regions %>% rename(!!name := tmp) 
  }
#Add the MDS coordinates to each window
  for (i in 1:40){
    j = i + 5
    win.regions[,j] <- mds.coords[,i]
  }
#Make plots of all the MDS to visualize patterns
pdf("chr7_landEU.MDSplots2.pdf", height=16,width=25)
print(
    win.regions %>%
      gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
      ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
      theme_bw()
  )
dev.off()
saveRDS(win.regions, file = "chr7_landEU.windows.rds")
#####这是做的画图的各种测试
pdf("chr7_landEU.MDSplots3.pdf", height=16,width=25)
win.regions %>%
  gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
  ggplot(aes(x=mid,y=value)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off()
pdf("chr7_landEU.MDSplots4.pdf", height=16,width=25)
win.regions %>%
  mutate(mds01_mean=rollapply(mds02,10,mean,align='right',fill=NA)) %>%
  ggplot(.,aes(x=mid,y=mds02)) + geom_point(alpha=0.5) +
  geom_line(aes(x=mid,y=mds01_mean))
dev.off()
pdf("chr7_landEU.MDSplots5.pdf", height=16,width=25)
win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-20)]) %>% 
  filter(mds == "mds03") %>%
  ggplot(aes(x=start,y=value)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off()
########现在是为了把所有的点找出来
win.fn.snp <- vcf_windower("chr7_landEU.withBarley.vcf.gz", size=100, type="snp", sites=sites) 
windows <- win.regions %>% filter(mds01 < -0.4) %>% pull(n) 
pca.test <- cov_pca(win.fn.snp(windows),k=2)  
win.regions %>% filter(n %in% windows) %>% summarize(start = min(start), end = max(end)) -> tmp.region
tmp.region$start
out <- pca.test
matrix.out <- t(matrix(out[4:length(out)],ncol=nrow(samples),byrow=T))
out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as.tibble() 
colnames(out) <- pull(samples)
###
out <- as_tibble(cbind(nms = names(out), t(out))) %>% 
  rename(name=nms,PC1=V2,PC2=V3) %>% 
  mutate(PC1 = as.double(PC1), PC2 = as.double(PC2))
#aaa = paste("bcftools query -H -f '%END [ %GT]\n' -r 8:", tmp.region$start,"-",tmp.region$end," ","chr7_landEU.withBarley.vcf.gz",'| sed s/\\#\\ //g |  sed s/\\:GT//g | sed s/END/pos/g > tmp.geno.txt',sep="")
system(paste("bcftools query -H -f '%END [ %GT]\n' -r 7:", tmp.region$start,"-",tmp.region$end," ","chr7_landEU.withBarley.vcf.gz",
             '| sed s/\\#\\ //g |  sed s/\\:GT//g | sed s/END/pos/g > tmp.geno.txt',sep=""))
##
read_delim("tmp.geno.txt",delim=" ",col_names = c("pos","blank", as.vector(samples$sample)),skip=1) %>%
  select(-blank) %>% mutate_if(.,  
                               is.character, 
                               str_replace_all, pattern = '0/0', replacement = "0") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '1/1', replacement = "2") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '0/1', replacement = "1") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = './.', replacement = "NA") -> snps
##现在做的是基本生物学统计量的分析
##这是杂合率
snps %>%  group_by(pos) %>%gather("name","genotype",1:(ncol(snps)-1)) %>%
  group_by(name, genotype) %>%
  summarize(count=n()) %>%
  spread(genotype, count) %>%
  summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity
pdf("chr7_landEU.MDSplots_het.pdf", height=8,width=8)
out %>%
  inner_join(.,heterozygosity) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=het)) + theme_bw() + scale_colour_viridis_c()
dev.off()
##这是maf
MAF <- 0.1
snps %>%  group_by(pos) %>%gather("name","genotype",2:(ncol(snps))) %>%
  group_by(pos,genotype) %>% summarize(count=n()) %>%
  spread(genotype, count,fill="0") %>% mutate(AA = as.numeric(`0`), Aa = as.numeric(`1`),aa = as.numeric(`2`)) %>%
  mutate(total = AA + Aa + aa, maf = ((aa * 2) + Aa)/(total*2)) %>%
  filter(maf > MAF, maf < (1-MAF)) %>% pull(pos) -> sites.maf
pdf("chr7_landEU.MDSplots_maf.pdf", height=8,width=8)
snps %>%  group_by(pos) %>%gather("name","genotype",2:(ncol(snps))) %>%
  filter(pos %in% sites.maf) %>% 
  inner_join(.,out) %>%
  filter(genotype != "NA") %>%
  ggplot(.,aes(x=as.character(pos),y=fct_reorder(name,PC1), fill=as.factor(genotype))) + geom_tile() +
  scale_fill_brewer(palette = "Set1")
dev.off()
##LD
#vcftools --gzvcf chr7_landEU.withBarley.vcf.gz --geno-r2 --ld-window-bp 500 --out ld_windows_500
#aaa = paste("bcftools query -H -f '%END [ %GT]\n' -r 2:", tmp.region$start,"-",tmp.region$end," ","chr7_landEU.withBarley.vcf.gz", "| vcftools --vcf - --stdout --maf 0.1 --geno-r2 > tmp.ld.txt",sep="")
#system(paste("bcftools query -H -f '%END [ %GT]\n' -r 7:", tmp.region$start,"-",tmp.region$end," ","chr7_landEU.withBarley.vcf.gz", " > tmp.ld1.vcf",sep=""))
system("vcftools --gzvcf chr7_landEU.withBarley.vcf.gz --geno-r2 --ld-window-bp 500000 --out tmp.ld.txt")
read_tsv("tmp.ld.txt.geno.ld") -> tmp.ld 
tmp.ld = na.omit(tmp.ld )
my_breaks <- c(0, 0.01, 0.1, 1)
pdf("chr7_landEU.MDSplots_ld.pdf", height=16,width=16)
tmp.ld %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=`R^2`)) + geom_tile() +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  scale_y_discrete(breaks = levels(as.factor(tmp.ld$POS1))[seq(1, length(levels(as.factor(tmp.ld$POS1))),50)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = levels(as.factor(tmp.ld$POS1))[seq(1, length(levels(as.factor(tmp.ld$POS1))),50)]) +
  xlab("position") + ylab("position")
dev.off()
################Trying to dplyr::select outliers for each MDS PC
#win.regions <- readRDS("chr7_landEU.windows.rds")
mds_pcs <- colnames(win.regions)[5:(ncol(win.regions)-1)]
min_windows <- 3
mds_clustering <- tibble(mds_coord = character(),direction = character(),clust_pvalue = numeric(), outliers = numeric(),n1_outliers=numeric(),
                         high_cutoff = numeric(),lower_cutoff=numeric(),chr=character())
for (mds_chosen in mds_pcs){
  print(paste("Processing",mds_chosen))
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    mutate(sd_mds = sd(the_mds)) %>%
    filter(the_mds > (sd_mds *3)) -> pos_windows
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    summarize(sd_mds = sd(the_mds)) %>% pull() -> sd_mds
  mds_high_cutoff <- sd_mds * 3
  mds_low_cutoff <- sd_mds * 2
  n_permutations <- 10
  if (nrow(pos_windows) >= min_windows){
    permutations <- matrix( nrow = n_permutations, ncol = 1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      #max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    
    x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow() 
    x <- x+1
    pvalue <- x/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("pos"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(mds_high_cutoff),lower_cutoff=as.numeric(mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("pos"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
    
  }
 
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    mutate(sd_mds = sd(the_mds)) %>%
    filter(the_mds < -(sd_mds *3)) -> neg_windows
  
  if (nrow(neg_windows) >= min_windows){
    permutations <- matrix( nrow = n_permutations, ncol = 1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(neg_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    
    x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow() 
    x <- x+1
    pvalue <- x/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(-mds_high_cutoff),lower_cutoff=as.numeric(-mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
    
  }
}
mds_clustering %>% filter(clust_pvalue < 2) -> sig_mds_clusters
#save(sig_mds_clusters,file="sig_mds_clusters.Rdata")
#load("sig_mds_clusters.Rdata")
outlier_windows <- tibble(chrom=character(),start=numeric(),end=numeric(),mid=numeric(),the_mds=numeric(),mds_coord=character(),outlier=character(),n=numeric())  
cluster_genotypes <- tibble(mds_coord=character(),name=character(),PC1=numeric(),genotype=character())
#For each mds outlier set that is chromosomally clustered, print a bunch of stuff about it.
for (i in 1:nrow(sig_mds_clusters)){
  coord <- pull(sig_mds_clusters[i,1])
  direction <- pull(sig_mds_clusters[i,2])
  high_cutoff <- pull(sig_mds_clusters[i,6])
  low_cutoff <- pull(sig_mds_clusters[i,7])
  cluster_chr <- pull(sig_mds_clusters[i,8])
  coord_direction <- paste(coord, "-",direction,sep="")
  print(paste("Testing",coord_direction))
    
  if (direction == "pos"){
    current_windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      mutate(outlier = case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Outlier",
                                 TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      dplyr::select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord = coord_direction) %>%
      mutate(ahead_n = n - lag(n),behind_n = abs(n - lead(n))) %>%
      mutate(min_dist = pmin(ahead_n, behind_n,na.rm=T)) %>%
      filter(min_dist < 100 ) %>%
      select(-ahead_n, -behind_n, -min_dist)
    
    windows <- current_windows %>% pull(n)
    
    outlier_windows <- rbind( current_windows, outlier_windows)
    
  }else{
    
    current_windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      mutate(outlier = case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Outlier",
                                 TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      dplyr::select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord = coord_direction) %>%
      mutate(ahead_n = n - lag(n),behind_n = abs(n - lead(n))) %>%
      mutate(min_dist = pmin(ahead_n, behind_n,na.rm=T)) %>%
      filter(min_dist < 100 ) %>%
      select(-ahead_n, -behind_n, -min_dist)
    
    windows <- current_windows %>% pull(n)
    
    outlier_windows <- rbind(current_windows, outlier_windows)  
  }
  out <- cov_pca(win.fn.snp(windows),k=2)
  matrix.out <- t(matrix(out[4:length(out)],ncol=nrow(samples),byrow=T))
  out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as.tibble() 
  colnames(out) <- pull(samples)
  
  out <- as_tibble(cbind(nms = names(out), t(out))) %>% 
    rename(name=nms,PC1=V2,PC2=V3) %>% 
    mutate(PC1 = as.double(PC1), PC2 = as.double(PC2))
  
  try_3_clusters <- try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]),(min(matrix.out[,1])+max(matrix.out[,1]))/2,max(matrix.out[,1]))))
  #try_3_clusters <- try(kmeans(matrix.out[,1], 3))
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]),max(matrix.out[,1])))
    out.normal <- out
    out.rotation <- out
  }else{
    rotation.ss <- tibble(rotation = numeric(),betweenss = numeric())
    for (i in seq(0.02, 3.14, 0.02)){
      rotated.matrix <- Rotation(matrix.out,i)
      
      rotated.kmeans <- kmeans(rotated.matrix[,1], 3, centers=c(min(rotated.matrix[,1]),(min(rotated.matrix[,1])+max(rotated.matrix[,1]))/2,max(rotated.matrix[,1])))
      rotated.tibble <- tibble(PC1 = as.numeric(rotated.matrix[,1]),
                               PC2 = as.numeric(rotated.matrix[,2]),
                               cluster = rotated.kmeans$cluster)
      
      # ggplot(rotated.tibble,aes(x=PC1,y=PC2,color=as.factor(cluster))) + 
      #  geom_point() + scale_color_brewer(palette = "Set1") + ggtitle(paste("BetweenSS =",round(rotated.kmeans$betweenss,3)))
      tmp.tibble <- tibble(rotation = as.numeric(i),betweenss = as.numeric(rotated.kmeans$betweenss))
      rotation.ss <- rbind(rotation.ss, tmp.tibble)
    }
    optimal.rotation <- rotation.ss[which(rotation.ss$betweenss == max(rotation.ss$betweenss)),1]$rotation
    out.normal <- out
    out.rotation <- out
    out.rotation[,2:3] <- Rotation(matrix.out,optimal.rotation)
    kmeans_cluster <-kmeans(out.rotation[,2], 3, centers=c(min(out.rotation[,2]),
                                                           (min(out.rotation[,2])+max(out.rotation[,2]))/2,
                                                           max(out.rotation[,2])))
    
  }
  
  out.rotation$cluster <- kmeans_cluster$cluster - 1 
  out.rotation$cluster <- as.character(out.rotation$cluster)
  out.rotation$mds_coord <- paste(coord, direction,sep="-")
  
  genotype.out <- out.rotation %>% dplyr::select(mds_coord,name,PC1,cluster) %>%
    rename(genotype = cluster) 
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
 
}

#################Way to visualize all outlier positions
#展示outlier
pdf("chr7_landEU.MDSplots_outlier_windows.pdf", height=16,width=16)
outlier_windows %>% 
  ggplot(.,aes(x=mid,y=mds_coord)) + geom_point() + 
  facet_wrap(~chrom,scales = "free_y")
dev.off()
#Clustering of outlier windows
mds_distances <- tibble(mds_coord = character(),mean_dist = numeric())
for (mds in unique(outlier_windows$mds_coord)){
  tmp_outliers <- outlier_windows %>%
    filter(mds_coord == mds)
  distances <- vector()
  for (i in 1:(nrow(tmp_outliers)-1)){
    for (j in (i+1):nrow(tmp_outliers)){
      dist <- abs(tmp_outliers$mid[i] - tmp_outliers$mid[j])
      distances <- c(distances, dist)
    }   
  }
  tmp_tibble <- tibble(mds_coord = as.character(mds),mean_dist = as.numeric(mean(distances)/1000000))
  mds_distances <- rbind(mds_distances, tmp_tibble)
}
#Counting of outlier windows
mds_counts <- tibble(mds_coord = character(),n_outliers = numeric())
for (mds in unique(outlier_windows$mds_coord)){
  count_outliers <- outlier_windows %>%
    filter(mds_coord == mds) %>% nrow()
  tmp_tibble <- tibble(mds_coord = as.character(mds),n_outliers = as.numeric(count_outliers))
  mds_counts <- rbind(mds_counts, tmp_tibble)
}
#Correlation between MDS for collapsing them
correlated_mds <- tibble(mds1 = character(),mds2 = character(),correlation = numeric())
for (mds1 in unique(outlier_windows$mds_coord)){
  for (mds2 in unique(outlier_windows$mds_coord)){
    if (mds1 == mds2){next}
    chr1 <- outlier_windows %>% filter(mds_coord == mds1) %>% dplyr::select(chrom) %>% unique() %>% pull()
    chr2 <- outlier_windows %>% filter(mds_coord == mds2) %>% dplyr::select(chrom) %>% unique() %>% pull()
    if (chr1 != chr2){next;}
    cluster_genotypes %>% mutate(mds_coord = gsub("_","-",mds_coord)) %>% filter( mds_coord == mds1 | mds_coord == mds2) %>%
      dplyr::select(-PC1) %>%
      spread(mds_coord, genotype) %>%  dplyr::select(-name) -> tmp
    x <- tmp %>% pull(1) %>% as.numeric()
    y <- tmp %>% pull(2) %>% as.numeric()
    test_result <- cor.test(x,y,na.rm=T)
    tmp_tibble <- tibble(mds1 = as.character(mds1),mds2 = as.character(mds2),correlation = as.numeric(abs(test_result$estimate)))
    correlated_mds <- rbind(correlated_mds, tmp_tibble)
    print(paste(mds2, test_result$estimate))
  }
}
#Check pairs with high correlation and pull out the mds that has fewer outlier windows.
total_mds_coords <- unique(outlier_windows$mds_coord)
min_cor <- 0.9
for (i in 1:nrow(correlated_mds)){
  if (correlated_mds[i,3] >= min_cor){
    count1 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,1])) %>% pull(n_outliers)
    count2 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,2])) %>% pull(n_outliers)
    if (count1 < count2){
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,1]))] -> total_mds_coords
    }else if(count1 >= count2){
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,2]))] -> total_mds_coords
    }
  }
}

##########################################################################################################################Redo plotting with non-redundent coordinates
outlier_windows <- tibble(chrom=character(),start=numeric(),end=numeric(),mid=numeric(),the_mds=numeric(),mds_coord=character(),outlier=character(),n=numeric())  
cluster_genotypes <- tibble(mds_coord=character(),name=character(),PC1=numeric(),genotype=character())
inversion_stats <- tibble(mds_coord=character(),betweenSS=numeric(),het_pvalue=numeric(),window_cluster=numeric(),
                          location=character())
#Prepare map data
  usa <- map_data('state')
  states <- map_data("state")
  #target_state <- map_data('state')
  target_state <- map_data('world')
  lat_range <- c(-25, 70)
  long_range <- c(-10,180)
  pie_size <- 2
	#labels <- read.table("test_land_lat.txt",header = T)
  pop_loc <- read.table("test_land_lat.txt",header = T)
  pop_loc %>% rename(population = name) %>% inner_join(.,labels) -> labels
##For each mds outlier set that is chromosomally clustered, print a bunch of stuff about it.
max_distance_between_outliers = 100
pdf("chr7_landEU.mdsoutliers.nonredundant.pdf",height=8,width=18)
for (i in 1:nrow(sig_mds_clusters)){
  coord <- pull(sig_mds_clusters[i,1])
  direction <- pull(sig_mds_clusters[i,2])
  if(! paste(coord,"-",direction,sep="") %in% total_mds_coords){
    print(paste("Skipping",coord, direction))
    next;
  }else{
    print(paste("Printing",coord, direction))
  }
  high_cutoff <- pull(sig_mds_clusters[i,6])
  low_cutoff <- pull(sig_mds_clusters[i,7])
  cluster_chr <- pull(sig_mds_clusters[i,8])
  coord_direction <- paste(coord, "-",direction,sep="")
  
  #Select outlier windows.
  if (direction == "pos"){
    current_windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      mutate(outlier = case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Outlier",
                                 TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      dplyr::select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord = coord_direction) %>%
      mutate(ahead_n = n - lag(n),behind_n = abs(n - lead(n))) %>%
      mutate(min_dist = pmin(ahead_n, behind_n,na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers ) %>%
      select(-ahead_n, -behind_n, -min_dist)
    
    windows <- current_windows %>% pull(n)
    
    outlier_windows <- rbind( current_windows, outlier_windows)
    
  }else{
    
    current_windows <- win.regions %>%
      mutate_(the_mds = coord ) %>% 
      mutate(outlier = case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Outlier",
                                 TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      dplyr::select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord = coord_direction) %>%
      mutate(ahead_n = n - lag(n),behind_n = abs(n - lead(n))) %>%
      mutate(min_dist = pmin(ahead_n, behind_n,na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers ) %>%
      select(-ahead_n, -behind_n, -min_dist)
    
    windows <- current_windows %>% pull(n)
    
    outlier_windows <- rbind(current_windows, outlier_windows)  
  }
  
  genome_plot <- win.regions %>%
    mutate_(the_mds = coord ) %>% 
    mutate(outlier = case_when(n %in% current_windows$n ~ "Outlier",
                               TRUE ~ "Non-outlier")) %>%
    ggplot(.,aes(x=mid/1000000,y=the_mds,color=outlier)) + geom_point() + facet_wrap(~chrom,scales= "free_x",nrow=1) +
    theme_bw() + scale_color_manual(values=c("black","#E41A1C")) + xlab("MB") + ylab(paste(coord))
  out <- cov_pca(win.fn.snp(windows),k=2)
  matrix.out <- t(matrix(out[4:length(out)],ncol=nrow(samples),byrow=T))
  out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as.tibble() 
  colnames(out) <- pull(samples)
  
  out <- as_tibble(cbind(nms = names(out), t(out))) %>% 
    rename(name=nms,PC1=V2,PC2=V3) %>% 
    mutate(PC1 = as.double(PC1), PC2 = as.double(PC2))
  
  try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]),(min(matrix.out[,1])+max(matrix.out[,1]))/2,max(matrix.out[,1]))))
  
  
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]),max(matrix.out[,1])))
    out.normal <- out
    out.rotation <- out
  }else{
    rotation.ss <- tibble(rotation = numeric(),betweenss = numeric())
    for (i in seq(0.02, 3.14, 0.02)){
      rotated.matrix <- Rotation(matrix.out,i)
      
      rotated.kmeans <- kmeans(rotated.matrix[,1], 3, centers=c(min(rotated.matrix[,1]),(min(rotated.matrix[,1])+max(rotated.matrix[,1]))/2,max(rotated.matrix[,1])))
      rotated.tibble <- tibble(PC1 = as.numeric(rotated.matrix[,1]),
                               PC2 = as.numeric(rotated.matrix[,2]),
                               cluster = rotated.kmeans$cluster)
      
      # ggplot(rotated.tibble,aes(x=PC1,y=PC2,color=as.factor(cluster))) + 
      #  geom_point() + scale_color_brewer(palette = "Set1") + ggtitle(paste("BetweenSS =",round(rotated.kmeans$betweenss,3)))
      tmp.tibble <- tibble(rotation = as.numeric(i),betweenss = as.numeric(rotated.kmeans$betweenss))
      rotation.ss <- rbind(rotation.ss, tmp.tibble)
    }
    optimal.rotation <- rotation.ss[which(rotation.ss$betweenss == max(rotation.ss$betweenss)),1]$rotation
    out.normal <- out
    out.rotation <- out
    out.rotation[,2:3] <- Rotation(matrix.out,optimal.rotation)
    kmeans_cluster <-kmeans(out.rotation[,2], 3, centers=c(min(out.rotation[,2]),
                                                           (min(out.rotation[,2])+max(out.rotation[,2]))/2,
                                                           max(out.rotation[,2])))   
  }
	pca_plot_normal <- out.normal %>%
    ggplot(.,aes(x=PC1,y=PC2)) + geom_point() + theme_bw()
  pca_plot_rotate <- out.rotation %>%
    ggplot(.,aes(x=PC1,y=PC2)) + geom_point() + theme_bw()
  
  out.rotation$cluster <- kmeans_cluster$cluster - 1 
  out.rotation$cluster <- as.character(out.rotation$cluster)
  out.rotation$mds_coord <- paste(coord, direction,sep="_")
  
  genotype.out <- out.rotation %>% dplyr::select(mds_coord,name,PC1,cluster) %>%
    rename(genotype = cluster) 
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
  
  hist_plot <- out.rotation %>%
    mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)) %>%
    ggplot(.,aes(PC1)) + geom_histogram(aes(fill=as.character(cluster))) +
    theme_bw() + scale_fill_brewer(palette = "Set1",name="Cluster")
  
  
  win.fn.snp(windows) %>% as.tibble() -> snps
  colnames(snps) <- pull(samples)
  
  snps %>% gather("name","genotype",1:ncol(snps)) %>%group_by(name, genotype) %>%
    summarize(count=n()) %>%
    spread(genotype, count) %>%
    summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity
  
  het_plot <- inner_join(out.rotation, heterozygosity) %>% 
    ggplot(.,aes(x=as.character(cluster),y=het,fill=as.character(cluster))) + 
    geom_boxplot() + scale_fill_brewer(palette = "Set1",name="Cluster") + theme_bw() + xlab("Cluster") + ylab("Heterozygosity")
  if(length(unique(kmeans_cluster$cluster)) == 3){
    map_plot <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=inner_join(out.rotation, labels) %>% 
                        group_by(population, lat, long, cluster) %>% 
                        tally() %>%
                        spread(., cluster, n,fill=0),
                      aes(x=long, y=lat, r=pie_size), 
                      cols=c("0","1","2"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0))
  }else{
    map_plot <- ggplot(target_state, aes(long, lat)) +
      geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
      coord_quickmap() +
      geom_scatterpie(data=inner_join(out.rotation, labels) %>% 
                        group_by(population, lat, long, cluster) %>% 
                        tally() %>%
                        spread(., cluster, n,fill=0),
                      aes(x=long, y=lat, r=2), 
                      cols=c("0","1"), color=NA, alpha=.8) +
      scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(limits = long_range, expand = c(0, 0)) +
      scale_y_continuous(limits = lat_range, expand = c(0, 0))
  }
  #Calculate stats about this potential inversion:
  #Window clustering.
  window.cluster <- current_windows %>%
    mutate(n.difference = n - lag(n)) %>% 
    summarize(mean_window_distance = mean(n.difference,na.rm=T)) %>%
    pull() %>% round(.,2)
  #Heterozygosity, is it higher in the middle cluster?
  
  if("try-error" %in% class(try_3_clusters)){
    high_pvalue <- "NA"
  }else{
    het.try1 <- try(t.test(inner_join(out.rotation, heterozygosity) %>% filter(cluster == 0) %>% pull(het) ,
                           inner_join(out.rotation, heterozygosity) %>% filter(cluster == 1) %>% pull(het)))
    het.try2 <- try(t.test(inner_join(out.rotation, heterozygosity) %>% filter(cluster == 2) %>% pull(het) ,
                           inner_join(out.rotation, heterozygosity) %>% filter(cluster == 1) %>% pull(het)))
    if(("try-error" %in% class(het.try1)) | ("try-error" %in% class(het.try2)) ){
      high_pvalue <- "NA"
    }else{
      het.test1 <- t.test(inner_join(out.rotation, heterozygosity) %>% filter(cluster == 0) %>% pull(het) ,
                          inner_join(out.rotation, heterozygosity) %>% filter(cluster == 1) %>% pull(het))
      het.test2 <- t.test(inner_join(out.rotation, heterozygosity) %>% filter(cluster == 2) %>% pull(het) ,
                          inner_join(out.rotation, heterozygosity) %>% filter(cluster == 1) %>% pull(het))
      
      if (het.test2$statistic < 0 & het.test1$statistic < 0){
        high_pvalue <- signif(max(het.test2$p.value,het.test1$p.value),3)
      }else{
        high_pvalue <- "NA"
      }
    }
  }
  #PC1 clustering
  PC1_cluster <- kmeans_cluster$betweenss
  
  #Find middle of inversion
  middle <- round(current_windows %>% summarize(middle_mid = median(mid)) %>% pull()/1000000)
  chromosome <- current_windows %>% head(1) %>% pull(chrom)
  
  #Save stats about each possible inversion
  tmp_stats <- tibble(mds_coord=as.character(coord_direction),betweenSS=as.numeric(PC1_cluster),
                      het_pvalue=as.numeric(high_pvalue),window_cluster=as.numeric(window.cluster),
                      location=as.character(paste(chromosome,":",middle,"MB", sep="")))
  inversion_stats <- rbind(inversion_stats, tmp_stats)
  print(
    grid.arrange(
      pca_plot_normal, pca_plot_rotate, hist_plot, het_plot ,map_plot, genome_plot,
      widths = c(1, 1, 1, 1, 1),
      layout_matrix = rbind(c(1, 2, 3, 4, 5),
                            c(6, 6, 6, 6, 6)),
      top = textGrob(paste(coord, "-", direction, " BetweenSS=",round(PC1_cluster,3),
                           " Het_pvalue=",high_pvalue," Window_cluster=",window.cluster,
                           " Location=",chromosome,":",middle,"MB", sep=""),gp=gpar(fontsize=20,font=1))
    )
  )
}
dev.off()

#Save cluster text output
write_tsv(outlier_windows, "chr7_landEU.mds_cluster_windows.txt")
write_tsv(cluster_genotypes, "chr7_landEU.mds_cluster_genotyped.txt")
write_tsv(inversion_stats,  "chr7_landEU.inversion_stats.txt")





