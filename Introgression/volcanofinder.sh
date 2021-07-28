#############################################volcano Finder#######################################
# VolcanoFinder is a program to perform genome-wide scans of adaptive introgression by using the
# spatial distribution of within-species polymorphism and between species substitution in the genome.
# input: 1.Allele frequency file 2.site frequency spectrum(sfs)3.ancestral file

###### prepare allele frequency file ######
# 1.choose target population of interest and parse vcf

# 2.simplify vcf to a table,
# ex:chr1	1145404	1145405	C	T	0.0261	0.0000	0.0856
zcat chr1.selected.vcf.gz |tail -n +22|awk '{split($8,a,";");print "chr1""\t"$2-1"\t"$2"\t"$4"\t"$5"\t"substr(a[7],5,length(a[7]))"\t"substr(a[8],9,length(a[8]))"\t"substr(a[9],8,length(a[9]))}' >chr1A_part1.bed

# 3.calculate allele numbers of ref and alt
# ex:chr1	821443	821444	T	1137	C	3
vcftools --gzvcf /data1/publicData/wheat/genotype/VMapII/VMapII/VMap2.0/chr001_vmap2.0.vcf.gz --counts --out chr1A_part1
cat chr1A_part1.frq.count |tail -n +2 |awk '{print "chr"$1"\t"$2-1"\t"$2"\t"substr($5,1,1)"\t"substr($5,3,length($5))"\t"substr($6,1,1)"\t"substr($6,3,length($6))}'>chr1A_part1.temp

# 4.extract and non-polarized sites to obtain Freq file
# ex:position x n folded
# 460000 9 100 0
# 460010 100 100 0
# 460210 30 78 1
bedtools intersect -a chr1.anc -b chr1A_part1.temp -wo |awk '{if($4==$8)print $2"\t"$11"\t"$11+$9"\t""0";else print $2"\t"$9"\t"$9+$11"\t""0" }' >polarized.freq
bedtools subtract -b chr1.anc -a chr1A_part1.temp |awk '{print $2"\t"$7"\t"$5+$7"\t""1"}' >>polarized.freq
# remove undistinguished sites
awk '{if($2==0)print $0 }' polarized.freq >x
awk '{if($2==$3&&$4==1)print $0 }' polarized.freq >>x
sort -k1,1n x>xx
awk 'FNR==NR{a[$0]++;next}(!($0 in a))' xx polarized.freq >FreqFile


# 5.sum to sfs Spect file
# 470 is alleles number
for var in {1..470};do awk -v j=$var '{if($2==j)print $2}' FreqFile|wc -l>>x ; done
# 471304005 is snp number
awk '{print NR"\t"$1/471304005}' x >SpectFile


###### run VolcanoFinder by blocks #######
# g:block length(ex:1000bp) G:step length(ex:1000snp)
#ex:VolcanoFinder -big 1000 FreqFile SpectFile -1 1 1 xbsj 1 50
./VolcanoFinder -big g FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK
./VolcanoFinder -bi G FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK

# merge results
#VolcanoFinder -m xbsj 50
./VolcanoFinder -m OutFile NBLOCK




#############################################volcano Finder#######################################
# VolcanoFinder is a program to perform genome-wide scans of adaptive introgression by using the
# spatial distribution of within-species polymorphism and between species substitution in the genome.
# input: 1.Allele frequency file 2.site frequency spectrum(sfs)3.ancestral file
##: /data1/home/xuebo/Projects/Speciation/volcanofinder

###### installation ######
wget -c http://degiorgiogroup.fau.edu/vf.html
tar -xzvf volcanofinder_v1.0.tar.gz
cd volcanofinder_v1.0
make
vim ~/.bashrc
source ~/.bashrc

###### prepare allele frequency file ######
# 1.calculate allele numbers of ref and al
# ex:chr1	821443	821444	T	1137	C	3
vcftools --gzvcf /data1/home/xuebo/Projects/Speciation/E3/chr36.all.vcf.gz  --counts --out chr36_part1
cat chr36_part1.frq.count |tail -n +2 |awk '{print $1"\t"$2-1"\t"$2"\t"substr($5,1,1)"\t"substr($5,3,length($5))"\t"substr($6,1,1)"\t"substr($6,3,length($6))}'>chr36_part1.temp
# 2.extract and non-polarized sites to obtain Freq file
# ex:position x n folded
# 460000 9 100 0
# 460010 100 100 0
# 460210 30 78 1
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' A_pos_ances.txt >  A_pos_ances.bed 
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' B_pos_ances.txt >  B_pos_ances.bed 
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' D_pos_ances.txt >  D_pos_ances.bed 
mkdir -p Dsplit_ances
for chr in `cut -f 1 D_pos_ances.bed | sort | uniq`; do
                grep -w $chr D_pos_ances.bed > Dsplit_ances/$chr.pos_ances.bed
done
***intersect,ancestralpos\tderived\tderived+ancestral
bedtools intersect -a ../anc/Dsplit_ances/36.pos_ances.bed  -b chr36_part1.temp -wo |awk '{if($4==$8)print $2"\t"$11"\t"$11+$9"\t""0";else print $2"\t"$9"\t"$9+$11"\t""0" }' >polarized.freq
***-apos
bedtools subtract -b ../anc/Dsplit_ances/36.pos_ances.bed  -a chr36_part1.temp |awk '{print $2"\t"$7"\t"$5+$7"\t""1"}' >>polarized.freq
# remove undistinguished sites
awk '{if($2==0)print $0 }' polarized.freq >x
awk '{if($2==$3&&$4==1)print $0 }' polarized.freq >>x
sort -k1,1n x>xx
awk 'FNR==NR{a[$0]++;next}(!($0 in a))' xx polarized.freq >FreqFile
sort -k1,1n FreqFile |sed '1i position\tx\tn\tfolded' > FreqFilesort
# 3.sum to sfs Spect file
# 604 is alleles number
for var in {1..604};do awk -v j=$var '{if($2==j)print $2}' FreqFilesort| wc -l >> tempp ; done
# 104072 is snp number  zcat /data1/home/xuebo/Projects/Speciation/E3/chr36.all.vcf.gz | grep -v "#" | wc -l 
# wc -l FreqFile
awk '{print NR"\t"$1/104072}' tempp >SpectFile

###### run VolcanoFinder by blocks #######
# g:block length(ex:1000bp) G:step length(ex:1000snp)
#ex:VolcanoFinder -big 1000 FreqFile SpectFile -1 1 1 xbsj 1 50
#VolcanoFinder -big g FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK
#VolcanoFinder -bi G FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK
VolcanoFinder -big 100 FreqFilesort SpectFile -1 1 1 xbsj 1 50

# merge results GSNP -big
#VolcanoFinder -m xbsj 50
VolcanoFinder -m OutFile NBLOCK
VolcanoFinder -bi 100 FreqFilesort SpectFile -1 1 1 test36.txt 1 50
VolcanoFinder -big 1000000 FreqFilesort SpectFile -1 1 1 xbsj 1 50


####1M
yafei@159.226.116.204:/data1/home/yafei/Project3/Finder/merge/find_1M
#
/data2/xuebo/Projects/Speciation/volcanofinder/merge1M
##GenWinsmooth
#!/usr/bin/Rscript
setwd("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M")
library(GenWin)
library(dplyr)
for(i in c(1:7)){
  file=paste("find_1M_chr",i,"A",sep="")
  Data1=read.table(file,sep="\t",header = T)
  Data2<-na.omit(Data1)
  Data <- filter(Data2, Data2[,2] != Inf)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  figure=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/pdf_",i,"A.pdf",sep="")
  pdf(figure,width=12,height=5)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,2]-mean(Data$LR))/sd(Data$LR)
  }
  NORM=splineAnalyze(Y=Z,map=Data$location,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
  dev.off()
  normScore=NORM$windowData
  outFile=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/soomth/normLR_chr",i,"A.txt",sep="")
  write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
}
for(i in c(1:7)){
  file=paste("find_1M_chr",i,"B",sep="")
  Data1=read.table(file,sep="\t",header = T)
  Data2<-na.omit(Data1)
  Data <- filter(Data2, Data2[,2] != Inf)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  figure=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/pdf_",i,"B.pdf",sep="")
  pdf(figure,width=12,height=5)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,2]-mean(Data$LR))/sd(Data$LR)
  }
  NORM=splineAnalyze(Y=Z,map=Data$location,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
  dev.off()
  normScore=NORM$windowData
  outFile=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/soomth/normLR_chr",i,"B.txt",sep="")
  write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
}
for(i in c(4:4)){
  file=paste("find_1M_chr","4","D",sep="")
  Data1=read.table(file,sep="\t",header = T)
  Data2<-na.omit(Data1)
  Data <- filter(Data2, Data2[,2] != Inf)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  figure=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/pdf_","4","D.pdf",sep="")
  pdf(figure,width=12,height=5)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,2]-mean(Data$LR))/sd(Data$LR)
  }
  NORM=splineAnalyze(Y=Z,map=Data$location,smoothness = 200,plotRaw=T,plotWindows = T,method = 4)
  dev.off()
  normScore=NORM$windowData
  outFile=paste("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/soomth/normLR_chr","4","D.txt",sep="")
  write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
}

###########D
sed 1d find_1M_chr1A | awk '$0="1A\t"$0' > addchr_1M_chr1A
sed 1d find_1M_chr2A | awk '$0="2A\t"$0' > addchr_1M_chr2A
sed 1d find_1M_chr3A | awk '$0="3A\t"$0' > addchr_1M_chr3A
sed 1d find_1M_chr4A | awk '$0="4A\t"$0' > addchr_1M_chr4A
sed 1d find_1M_chr5A | awk '$0="5A\t"$0' > addchr_1M_chr5A
sed 1d find_1M_chr6A | awk '$0="6A\t"$0' > addchr_1M_chr6A
sed 1d find_1M_chr7A | awk '$0="7A\t"$0' > addchr_1M_chr7A
sed 1d find_1M_chr1B | awk '$0="1B\t"$0' > addchr_1M_chr1B
sed 1d find_1M_chr2B | awk '$0="2B\t"$0' > addchr_1M_chr2B
sed 1d find_1M_chr3B | awk '$0="3B\t"$0' > addchr_1M_chr3B
sed 1d find_1M_chr4B | awk '$0="4B\t"$0' > addchr_1M_chr4B
sed 1d find_1M_chr5B | awk '$0="5B\t"$0' > addchr_1M_chr5B
sed 1d find_1M_chr6B | awk '$0="6B\t"$0' > addchr_1M_chr6B
sed 1d find_1M_chr7B | awk '$0="7B\t"$0' > addchr_1M_chr7B
sed 1d find_1M_chr1D | awk '$0="1D\t"$0' > addchr_1M_chr1D
sed 1d find_1M_chr2D | awk '$0="2D\t"$0' > addchr_1M_chr2D
sed 1d find_1M_chr3D | awk '$0="3D\t"$0' > addchr_1M_chr3D
sed 1d find_1M_chr4D | awk '$0="4D\t"$0' > addchr_1M_chr4D
sed 1d find_1M_chr5D | awk '$0="5D\t"$0' > addchr_1M_chr5D
sed 1d find_1M_chr6D | awk '$0="6D\t"$0' > addchr_1M_chr6D
sed 1d find_1M_chr7D | awk '$0="7D\t"$0' > addchr_1M_chr7D
cat add* > chr_1M_chrall.txt
awk '{print NR"\t"$0}' chr_1M_chrall.txt > chr_1M_chrall_addline.txt 
##R
  Data=read.table("chr_1M_chrall_addline2.txt",sep="\t",header = F)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  #pdf("/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/pdf_all",width=12,height=5)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,4]-mean(Data$V4))/sd(Data$V4)
  }
  NORM=splineAnalyze(Y=Z,map=Data$V1,smoothness = 2,plotRaw=T,plotWindows = T,method = 4)
  #dev.off()
  normScore=NORM$windowData
  write.table(normScore,"/data2/xuebo/Projects/Speciation/volcanofinder/merge1M/soomth/normLR_chrall.txt",sep="\t",col.names = T,row.names = F)








###haplotype block
#https://people.maths.bris.ac.uk/~madjl/finestructure/toolsummary.html
wget -c https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter-0.0.4.tar.gz
tar -xzvf chromopainter-0.0.4.tar.gz
cd chromopainter-0.0.4/
./configure
make
chromopainter -help
wget -c https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromocombine-0.0.4.tar.gz
tar -xzvf chromocombine-0.0.4.tar.gz 
cd chromocombine-0.0.4
./configure
make


##
./ChromoPainter -g [haplotype inle] -r [recom rate inle] -f [donor list inle] -o [output lename]
#input
(1) the SNP data for a set of admixed recipient chromosomes, 
(2) the SNP data for a set of donor chromosomes thought to represent the sources of admixture in the recipient chromosomes, 
(3) a genetic-map representing the recombination distance between each pair of contiguous SNPs. #-u

##
Q gene
196896702-196900381 
Q5K5K    196891702-196905381
#
bcftools view -r 26:196896702-196900381 -Oz /data1/home/xuebo/Projects/Speciation/E3/chr26.all.vcf.gz -o chr26Q.vcf.gz &
java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=chr26Q.vcf.gz out=chr26Q.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001 
#vcfplink format
vcftools --gzvcf chr26Q.beagle.vcf.gz --plink --out chr26Q.beagle 
plink --file chr26Q.beagle --recode12 --out chr26Q.beagle.plink
cut -d " " -f 1 chr26Q.beagle.plink.ped > name.txt
./perlcode/plink2chromopainter.pl -p=chr26Q.beagle.plink.ped -m=chr26Q.beagle.plink.map -o=chr26Q.beagle.plink.phase -d=name.txt
##
./perlcode/makeuniformrecfile.pl -in -i 10 chr26Q.beagle.plink.phase chr26Q.beagle.plink.inrec

chromopainter -g chr26Q.beagle.plink.phase -r chr26Q.beagle.plink.inrec -o chr26Q.cp

chromopainter -g chr26Q.beagle.plink.phase -r chr26Q.beagle.plink.inrec -J -f donorlist.txt -o chr26Q.cp
chromopainter -g chr26Q.beagle.plink.phase -r chr26Q.beagle.plink.inrec -J -f donorlist2.txt -o chr26Q.cp2




##top30block
java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=allA.vcf.gz out=allA.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001	&
java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=allB.vcf.gz out=allB.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001  &
java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=allD.vcf.gz out=allD.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001	&
vcftools --gzvcf allA.beagle.vcf.gz --plink --out allA.beagle
vcftools --gzvcf allB.beagle.vcf.gz --plink --out allB.beagle &
vcftools --gzvcf allD.beagle.vcf.gz --plink --out allD.beagle &
plink --file allA.beagle --recode12 --chr-set 42 --out allA.beagle.plink
plink --file allB.beagle --recode12 --chr-set 42 --out allB.beagle.plink
plink --file allD.beagle --recode12 --chr-set 42 --out allD.beagle.plink
##
cut -d " " -f 1 allA.beagle.plink.ped > nameA.txt
../perlcode/plink2chromopainter.pl -p=allA.beagle.plink.ped -m=allA.beagle.plink.map -o=allA.beagle.plink.phase -d=nameA.txt
cut -d " " -f 1 allB.beagle.plink.ped > nameB.txt
../perlcode/plink2chromopainter.pl -p=allB.beagle.plink.ped -m=allB.beagle.plink.map -o=allB.beagle.plink.phase -d=nameB.txt
cut -d " " -f 1 allD.beagle.plink.ped > nameD.txt
../perlcode/plink2chromopainter.pl -p=allD.beagle.plink.ped -m=allD.beagle.plink.map -o=allD.beagle.plink.phase -d=nameD.txt
##
../perlcode/makeuniformrecfile.pl allA.beagle.plink.phase allA.beagle.plink.inrec
../perlcode/makeuniformrecfile.pl allB.beagle.plink.phase allB.beagle.plink.inrec
../perlcode/makeuniformrecfile.pl allD.beagle.plink.phase allD.beagle.plink.inrec
##
nohup chromopainter -g allA.beagle.plink.phase -r allA.beagle.plink.inrec -in -i 2 -a 0 0  -f donorlist_A.txt -o allA.cpainter > nohupoutA 2>& 1 &
nohup chromopainter -g allB.beagle.plink.phase -r allB.beagle.plink.inrec -in -i 2 -a 0 0  -f donorlist_B.txt -o allB.cpainter > nohupoutB 2>& 1 &
nohup chromopainter -g allD.beagle.plink.phase -r allD.beagle.plink.inrec -in -i 2 -a 0 0  -f donorlist_D.txt -o allD.cpainter > nohupoutD 2>& 1 &

#############
*****/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockA_order
#!/bin/bash
a=1
cat /data1/home/xuebo/Projects/Speciation/volcanofinder/EUland/allA.bed | while read line
do
    echo $line > test111.txt     
    cat test111.txt | cut -d ":" -f 1 > nameChr.txt
		nameChr=$(sed -n 1p nameChr.txt)
    #echo $nameChr
    echo -e "bcftools view -r $line -Oz /data1/home/xuebo/Projects/Speciation/E3/chr${nameChr}.all.vcf.gz -o top30_blockA_${a}.vcf.gz &"  >> runbcftools.sh
    a=$(($a+1))
    echo $a
done
sh runbcftools.sh
#!/bin/bash
for i in {1..15}
do 
	vcftools --gzvcf top30_blockA_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderA.txt \
   --maf 0.0001  --recode --stdout | bgzip -c > top30_blockA_order_${i}.vcf.gz &
done
#!/bin/bash
for i in {1..15}
do 
	WGS --model vcf --type keep --file top30_blockA_order_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderA.txt \
	--out top30_blockA_order_${i}_2.vcf 
	bgzip top30_blockA_order_${i}_2.vcf 
done
#!/bin/bash
for i in {1..15}
do 
	java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=top30_blockA_order_${i}_2.vcf.gz out=top30_blockA_${i}.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001	
	vcftools --gzvcf top30_blockA_${i}.beagle.vcf.gz --plink --out top30_blockA_${i}.beagle 
	plink --file top30_blockA_${i}.beagle  --recode12 --chr-set 42 --out top30_blockA_${i}.beagle.plink
done
cut -d " " -f 1 top30_blockA_1.beagle.plink.ped > blockA.txt
#!/bin/bash
for i in {1..15}
do 
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/plink2chromopainter.pl -p=top30_blockA_${i}.beagle.plink.ped -m=top30_blockA_${i}.beagle.plink.map -o=top30_blockA_${i}.beagle.plink.phase -d=blockA.txt
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/makeuniformrecfile.pl top30_blockA_${i}.beagle.plink.phase top30_blockA_${i}.beagle.plink.inrec
done
#!/bin/bash
for i in {1..15}
do 
	chromopainter -g top30_blockA_${i}.beagle.plink.phase -r top30_blockA_${i}.beagle.plink.inrec -in -i 1 -a 0 0  -f /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/donor_orderA.txt -o blockA_${i}.cpainter  &
done
nohup sh getcpainter.sh > nohupout 2>& 1 &

*****/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockB_order
#!/bin/bash
a=1
cat /data1/home/xuebo/Projects/Speciation/volcanofinder/EUland/allB.bed | while read line
do
    echo $line > test111.txt     
    cat test111.txt | cut -d ":" -f 1 > nameChr.txt
		nameChr=$(sed -n 1p nameChr.txt)
    #echo $nameChr
    echo -e "bcftools view -r $line -Oz /data1/home/xuebo/Projects/Speciation/E3/chr${nameChr}.all.vcf.gz -o top30_blockB_${a}.vcf.gz &"  >> runbcftools.sh
    a=$(($a+1))
    echo $a
done
sh runbcftools.sh
#!/bin/bash
for i in {1..8}
do 
	vcftools --gzvcf top30_blockB_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderB.txt \
   --maf 0.0001  --recode --stdout | bgzip -c > top30_blockB_order_${i}.vcf.gz &
done
#!/bin/bash
for i in {1..8}
do 
	WGS --model vcf --type keep --file top30_blockB_order_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderB.txt \
	--out top30_blockB_order_${i}_2.vcf 
	bgzip top30_blockB_order_${i}_2.vcf 
done
#!/bin/bash
for i in {1..8}
do 
	java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=top30_blockB_order_${i}_2.vcf.gz out=top30_blockB_${i}.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001	
	vcftools --gzvcf top30_blockB_${i}.beagle.vcf.gz --plink --out top30_blockB_${i}.beagle 
	plink --file top30_blockB_${i}.beagle  --recode12 --chr-set 42 --out top30_blockB_${i}.beagle.plink
done
cut -d " " -f 1 top30_blockB_1.beagle.plink.ped > blockB.txt
#!/bin/bash
for i in {1..8}
do 
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/plink2chromopainter.pl -p=top30_blockB_${i}.beagle.plink.ped -m=top30_blockB_${i}.beagle.plink.map -o=top30_blockB_${i}.beagle.plink.phase -d=blockB.txt
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/makeuniformrecfile.pl top30_blockB_${i}.beagle.plink.phase top30_blockB_${i}.beagle.plink.inrec
done
#!/bin/bash
for i in {1..8}
do 
	chromopainter -g top30_blockB_${i}.beagle.plink.phase -r top30_blockB_${i}.beagle.plink.inrec -in -i 1 -a 0 0  -f /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/donor_orderB.txt -o blockB_${i}.cpainter  &
done
nohup sh getcpainter.sh > nohupout 2>& 1 &

*****/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockD_order
#!/bin/bash
a=1
cat /data1/home/xuebo/Projects/Speciation/volcanofinder/EUland/allD.bed | while read line
do
    echo $line > test111.txt     
    cat test111.txt | cut -d ":" -f 1 > nameChr.txt
		nameChr=$(sed -n 1p nameChr.txt)
    #echo $nameChr
    echo -e "bcftools view -r $line -Oz /data1/home/xuebo/Projects/Speciation/E3/chr${nameChr}.all.vcf.gz -o top30_blockD_${a}.vcf.gz &"  >> runbcftools.sh
    a=$(($a+1))
    echo $a
done
sh runbcftools.sh
#!/bin/bash
for i in {1..6}
do 
	vcftools --gzvcf top30_blockD_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderD.txt \
   --maf 0.0001  --recode --stdout | bgzip -c > top30_blockD_order_${i}.vcf.gz &
done
#!/bin/bash
for i in {1..6}
do 
	WGS --model vcf --type keep --file top30_blockD_order_${i}.vcf.gz --keep /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/Cpainter_orderD.txt \
	--out top30_blockD_order_${i}_2.vcf 
	bgzip top30_blockD_order_${i}_2.vcf 
done
#!/bin/bash
for i in {1..6}
do 
	java -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=48 gt=top30_blockD_order_${i}_2.vcf.gz out=top30_blockD_${i}.beagle phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001	
	vcftools --gzvcf top30_blockD_${i}.beagle.vcf.gz --plink --out top30_blockD_${i}.beagle 
	plink --file top30_blockD_${i}.beagle  --recode12 --chr-set 42 --out top30_blockD_${i}.beagle.plink
done
cut -d " " -f 1 top30_blockD_1.beagle.plink.ped > blockD.txt
#!/bin/bash
for i in {1..6}
do 
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/plink2chromopainter.pl -p=top30_blockD_${i}.beagle.plink.ped -m=top30_blockD_${i}.beagle.plink.map -o=top30_blockD_${i}.beagle.plink.phase -d=blockD.txt
	/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/perlcode/makeuniformrecfile.pl top30_blockD_${i}.beagle.plink.phase top30_blockD_${i}.beagle.plink.inrec
done
#!/bin/bash
for i in {1..6}
do 
	chromopainter -g top30_blockD_${i}.beagle.plink.phase -r top30_blockD_${i}.beagle.plink.inrec -in -i 1 -a 0 0  -f /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/groupOrder/donor_orderD.txt -o blockD_${i}.cpainter  &
done
nohup sh getcpainter.sh > nohupout 2>& 1 &
*****ChromoCombine  
#/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockD_order
chromocombine -C -u -o ChromoCombine_top30_blockD -d /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockD_order/ChromoCombine
#/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockB_order
chromocombine -C -u -o ChromoCombine_top30_blockB -d /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockB_order/ChromoCombine
#/data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockA_order
chromocombine -C -u -o ChromoCombine_top30_blockA -d /data1/home/xuebo/Projects/Speciation/haplotype/chromopainter/top30_blockA_order/ChromoCombine



#####urartutop/data1/home/xuebo/Projects/Speciation/volcanofinder/GO_analysis
ls -l |awk '{print $NF}' > listA
#!/bin/bash
cat ../Agene/listA | while read line
do
	echo $line > test111.txt
	cat test111.txt | cut -d "." -f 1 > name.txt
	name=$(sed -n 1p name.txt)
	java -jar /data1/home/xuebo/Projects/Speciation/javaCode/GO_analysis.jar --file1 ../Agene/$line --file2 /data1/home/xuebo/Projects/Speciation/javaCode/mart_export.txt --out A_${name}_gene_GO.txt &
done
ls -l |awk '{print $NF}' > listB
#!/bin/bash
cat ../Bgene/listB | while read line
do
	echo $line > test111.txt
	cat test111.txt | cut -d "." -f 1 > name.txt
	name=$(sed -n 1p name.txt)
	java -jar /data1/home/xuebo/Projects/Speciation/javaCode/GO_analysis.jar --file1 ../Bgene/$line --file2 /data1/home/xuebo/Projects/Speciation/javaCode/mart_export.txt --out B_${name}_gene_GO.txt &
done
ls -l |awk '{print $NF}' > listD
#!/bin/bash
cat ../Dgene/listD | while read line
do
	echo $line > test111.txt
	cat test111.txt | cut -d "." -f 1 > name.txt
	name=$(sed -n 1p name.txt)
	java -jar /data1/home/xuebo/Projects/Speciation/javaCode/GO_analysis.jar --file1 ../Dgene/$line --file2 /data1/home/xuebo/Projects/Speciation/javaCode/mart_export.txt --out D_${name}_gene_GO.txt &
done

sort GO_B.txt | uniq -c | sort -k1,1n > Go_B_sort.txt


















