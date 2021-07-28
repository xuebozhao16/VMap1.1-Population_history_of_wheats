library(tidyverse)
library(lostruct)
library(colorspace)
library(RColorBrewer)
library(ggmap)
library(zoo)
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/lostruct")
samples <- read.table("sub_Urartu.txt",header =F)
colnames(samples) <- c("sample")
snps <- read_vcf("test_urartu.chr2.vcf.gz") #这是在跑bcftools
## save(snps,file="snps.Rdata")
## load("snps.Rdata")
eigenstuff <- eigen_windows(snps, win=100, k=2)  #这一步是多线程
windist <- pc_dist( eigenstuff, npc=2 )
windist_na <- which(is.na(windist), TRUE)
na.inds <- is.na(windist[,1]) 
if (sum(na.inds) == length(na.inds)){
  na.inds <- is.na(windist[,2]) 
}
mds <- cmdscale( windist[!na.inds,!na.inds], eig=TRUE, k=20 )
plot( mds$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(windist)) )
dev.off()
sites <- vcf_positions("test_urartu.chr2.vcf")
win.fn.snp <- vcf_windower("test_urartu.chr2.vcf", size=100, type="snp", sites=sites) 
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions$n <- 1:nrow(win.regions)
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
for (k in 1:40){
  str_pad(k, 2, pad = "0")
  
  name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions %>% rename(!!name := tmp) 
}
for (i in 1:20){
  j = i + 5
  win.regions[,j] <- mds.coords[,i]
}
pdf("test_chr2.MDSplots2.pdf", height=16,width=25)
print(
  win.regions %>%
    gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-20)]) %>% 
    ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
    theme_bw()
)
dev.off()
saveRDS(win.regions, file = "test_chr2.windows.rds")
#####这是做的画图的各种测试
pdf("test_chr2.MDSplots3.pdf", height=16,width=25)
win.regions %>%
  gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
  ggplot(aes(x=mid,y=value)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off()
pdf("test_chr2.MDSplots4.pdf", height=16,width=25)
win.regions %>%
  mutate(mds01_mean=rollapply(mds02,10,mean,align='right',fill=NA)) %>%
  ggplot(.,aes(x=mid,y=mds02)) + geom_point(alpha=0.5) +
  geom_line(aes(x=mid,y=mds01_mean))
dev.off()
pdf("test_chr2.MDSplots5.pdf", height=16,width=25)
win.regions %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-20)]) %>% 
  filter(mds == "mds03") %>%
  ggplot(aes(x=start,y=value)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off()
########
win.fn.snp <- vcf_windower("test_urartu.chr2.vcf.gz", size=100, type="snp", sites=sites) 
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
aaa = paste("bcftools query -H -f '%END [ %GT]\n' -r 2:", tmp.region$start,"-",tmp.region$end," ","test_urartu.chr2.vcf.gz",'| sed s/\\#\\ //g |  sed s/\\:GT//g | sed s/END/pos/g > tmp.geno.txt',sep="")
system(paste("bcftools query -H -f '%END [ %GT]\n' -r 2:", tmp.region$start,"-",tmp.region$end," ","test_urartu.chr2.vcf.gz",
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
##
snps %>%  group_by(pos) %>%gather("name","genotype",1:(ncol(snps)-1)) %>%
  group_by(name, genotype) %>%
  summarize(count=n()) %>%
  spread(genotype, count) %>%
  summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity
##
pdf("test_chr2.MDSplots_het.pdf", height=8,width=8)
out %>%
  inner_join(.,heterozygosity) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=het)) + theme_bw() + scale_colour_viridis_c()
dev.off()
##
MAF <- 0.1
snps %>%  group_by(pos) %>%gather("name","genotype",2:(ncol(snps))) %>%
  group_by(pos,genotype) %>% summarize(count=n()) %>%
  spread(genotype, count,fill="0") %>% mutate(AA = as.numeric(`0`), Aa = as.numeric(`1`),aa = as.numeric(`2`)) %>%
  mutate(total = AA + Aa + aa, maf = ((aa * 2) + Aa)/(total*2)) %>%
  filter(maf > MAF, maf < (1-MAF)) %>% pull(pos) -> sites.maf
pdf("test_chr2.MDSplots_maf.pdf", height=8,width=8)
snps %>%  group_by(pos) %>%gather("name","genotype",2:(ncol(snps))) %>%
  filter(pos %in% sites.maf) %>% 
  inner_join(.,out) %>%
  filter(genotype != "NA") %>%
  ggplot(.,aes(x=as.character(pos),y=fct_reorder(name,PC1), fill=as.factor(genotype))) + geom_tile() +
  scale_fill_brewer(palette = "Set1")
dev.off()
##
aaa = paste("bcftools query -H -f '%END [ %GT]\n' -r 2:", tmp.region$start,"-",tmp.region$end," ","test_urartu.chr2.vcf.gz", "| vcftools --vcf - --stdout --maf 0.1 --geno-r2 > tmp.ld.txt",sep="")
system(paste("bcftools query -H -f '%END [ %GT]\n' -r 2:", tmp.region$start,"-",tmp.region$end," ","test_urartu.chr2.vcf", "| vcftools --vcf - --stdout --maf 0.1 --geno-r2 > tmp.ld.txt",sep=""))
read_tsv("tmp.ld.txt") -> tmp.ld 
pdf("test_chr2.MDSplots_het.maf", height=8,width=8)
tmp.ld %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=`R^2`)) + geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  scale_y_discrete(breaks = levels(as.factor(tmp.ld$POS1))[seq(1, length(levels(as.factor(tmp.ld$POS1))),50)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = levels(as.factor(tmp.ld$POS1))[seq(1, length(levels(as.factor(tmp.ld$POS1))),50)]) +
  xlab("position") + ylab("position")
dev.off()
################Trying to dplyr::select outliers for each MDS PC
mds_pcs <- colnames(win.regions)[5:(ncol(win.regions)-1)]
min_windows <- 4
mds_clustering <- tibble(mds_coord = character(),direction = character(),clust_pvalue = numeric(), outliers = numeric(),n1_outliers=numeric(),
                         high_cutoff = numeric(),lower_cutoff=numeric(),chr=character())
for (mds_chosen in mds_pcs){
  print(paste("Processing",mds_chosen))
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    mutate(sd_mds = sd(the_mds)) %>%
    filter(the_mds > (sd_mds *4)) -> pos_windows
  win.regions %>%
    mutate_(the_mds = mds_chosen ) %>%
    summarize(sd_mds = sd(the_mds)) %>% pull() -> sd_mds
  mds_high_cutoff <- sd_mds * 4
  mds_low_cutoff <- sd_mds * 3
  n_permutations <- 1000
  if (nrow(pos_windows) >= min_windows){
    permutations <- matrix( nrow = n_permutations, ncol = 1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
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
    filter(the_mds < -(sd_mds *4)) -> neg_windows
  
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
mds_clustering %>% filter(clust_pvalue < 0.5) -> sig_mds_clusters
outlier_windows <- tibble(chrom=character(),start=numeric(),end=numeric(),mid=numeric(),the_mds=numeric(),mds_coord=character(),outlier=character(),n=numeric())  
cluster_genotypes <- tibble(mds_coord=character(),name=character(),PC1=numeric(),genotype=character())






