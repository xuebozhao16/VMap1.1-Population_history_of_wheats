##https://github.com/simonhmartin/genomics_general
####周姚的代码
1. data converting (vcf to geno)
WGS --model vcf --type toDstatistic --file ./data/genotype/A.vcf --out ./data/genotype/A.geno
bgzip A.geno
2. Generating scripts
java -jar ./code/java/WGS.jar --model GenerateScripts --type slicedD --file A.geno.gz --groupFile ./data/group/A.info.txt --popFile ../../group/all_Dstat.txt
3. run
python ../genomics_general/ABBABABAwindows.py  --windType sites -f phased -g Amaf.geno.gz -P1 sphaerococcum -P2 cultivar  -P3 urartu -O outgroup  \
-o Amaf9.csv.gz -w 100 --overlap 50 -T 1 -m 3 --popsFile ../../group/all_Dstat3.txt --writeFailedWindows &

python ../genomics_general/ABBABABAwindows.py  --windType sites -f phased -g AB.geno.gz -P1 domesticated_emmer -P2 bread_wheat  -P3 wild_emmer -O outgroup \
 -o AB1.csv.gz -w 100 --overlap 50 -T 1 -m 3 --popsFile ../../group/all_Dstat_081919.txt --writeFailedWindows &

###Anaconda解决Python2和python3 共存问题,之前204装的是3，但是这个软件要依赖2.7
https://foofish.net/compatible-py2-and-py3.html  学习网址
***基于 python3.6 创建一个名为test_py3 的环境
conda create --name test_py3 python=3.6 
***基于 python2.7 创建一个名为test_py2 的环境
conda create --name test_py2 python=2.7
***激活 test 环境
activate test_py2  # windows
source activate test_py2 # linux/mac
***切换到python3
activate test_py3

##因为至少加上两个outgroup,AB lineage使用tauschii和barley作为外类群，tauschii要做scan，tauschii往AB基因组上做的mapping
/data1/home/yaozhou/data/ref/wheat/genome/outgroup/bwa 这是周姚的数据存放的位置
/data2/xuebo/Projects/Speciation/introgression/outgroup/Ae  这是204上面的位置
cp.sh  
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do 
	cp /data2/xuebo/Projects/Speciation/hapScan/hapPos_ABD/chr${chr}.pos.txt  /data2/xuebo/Projects/Speciation/hapScan/hapPos/
	cp /data2/xuebo/Projects/Speciation/hapScan/posAllele_ABD/chr${chr}.allele.txt  /data2/xuebo/Projects/Speciation/hapScan/posAllele/
done
#之前文件里有一个空格
noblank.sh
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do 
	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/cutALTblank2.jar  --file1 /data2/xuebo/Projects/Speciation/tree/withBarley/row_chr${chr}.withBarley.vcf.gz \
	--out /data2/xuebo/Projects/Speciation/tree/withBarley/AD_noblank/row_chr${chr}.noblank_withBarley.vcf  &
	monitor java 20 60s
done
nohup sh noblank.sh > nohupnoblank 2>& 1 &
#################################################scan
***** runningAe.sh 这个是AB的外类群，所以只做AB
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
    java -jar -Xms10g -Xmx20g /data2/xuebo/Projects/Speciation/hapScan/hapScanner2_jar/PlantGenetics.jar \
    /data2/xuebo/Projects/Speciation/introgression/scanAe/parameters_file/parameters_hapScannerAe_chr${chr}.txt > logchr${chr}.txt 
done
nohup sh runningAe.sh > nohupAe 2>& 1 &
***** 压缩
#!/bin/bash
for chr in {"001","002","007","008","013","014","019","020","025","026","031","032","037","038","003","004","009","010","015","016","021","022","027","028","033","034","039","040"}
do
	  i=$(echo $chr | sed 's/^0*//g')
	  #echo $i
	  bgzip -c /data2/xuebo/Projects/Speciation/introgression/scanAe/outchr${i}/VCF/chr${chr}.vcf > /data2/xuebo/Projects/Speciation/introgression/scanAe/scan_Ae_vcf/chr${i}.Ae.vcf.gz &
done
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
	tabix chr${chr}.Ae.vcf.gz &
done
****合并，查看行数  行数正确
grep -v "#" /data2/xuebo/Projects/Speciation/introgression/scanAe/outchr1/VCF/chr001.vcf | wc -l  
grep -v "#" /data2/xuebo/Projects/Speciation/introgression/scanAe/outchr10/VCF/chr010.vcf | wc -l   
*merge.sh
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
		bcftools merge  -Oz /data2/xuebo/Projects/Speciation/tree/withBarley/row_chr${chr}.withBarley.vcf.gz  /data2/xuebo/Projects/Speciation/introgression/scanAe/scan_Ae_vcf/chr${chr}.Ae.vcf.gz \
    -o /data2/xuebo/Projects/Speciation/introgression/scanAe/raw_withbarleyAe/row_chr${chr}.withBarleyAe.vcf.gz &
done
nohup sh merge.sh > nohupmerge 2>& 1 &
*merge2.sh
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
		bcftools merge  -Oz /data2/xuebo/Projects/Speciation/tree/withBarley/AD_noblank/row_chr${chr}.noblank_withBarley.vcf.gz /data2/xuebo/Projects/Speciation/introgression/scanAe/scan_Ae_vcf/chr${chr}.Ae.vcf.gz \
    -o /data2/xuebo/Projects/Speciation/introgression/scanAe/raw_withbarleyAe/row_chr${chr}.withBarleyAe.vcf.gz &
done
nohup sh merge2.sh > nohupmerge2 2>& 1 &
*merge3.sh
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/C27_mergeTwoVCFbycol.jar --file1 /data2/xuebo/Projects/Speciation/introgression/scanAe/scan_Ae_vcf/chr${chr}.Ae.vcf.gz \
	--file2 /data2/xuebo/Projects/Speciation/tree/withBarley/AD_noblank/row_chr${chr}.noblank_withBarley.vcf --out /data2/xuebo/Projects/Speciation/introgression/scanAe/raw_withbarleyAe/row_chr${chr}.withBarleyAe.vcf &
done
nohup sh merge3.sh > nohupmerge3 2>& 1 &
zcat row_chr1.withBarleyAe.vcf | grep -v "#" | wc -l
zcat row_chr3.withBarleyAe.vcf.gz | grep -v "#" | wc -l
*tbi.sh
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
  tabix row_chr${chr}.withBarleyAe.vcf.gz &
done
###把VCF文件里面Ae和barley是./.的去掉，保留Ae和barley都是有数值的
getwithBarleyAe.sh 
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/C26_getwithBarley_Ae_VCF.jar --file1 /data2/xuebo/Projects/Speciation/introgression/scanAe/raw_withbarleyAe/row_chr${chr}.withBarleyAe.vcf.gz \
	--out /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr${chr}.withBarleyAe.vcf &
	monitor java 20 60s
done
nohup sh getwithBarleyAe.sh > nohup.out 2>& 1 &
getwithBarleyAe2.sh 
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/C26_getwithBarley_Ae_VCF.jar --file1 /data2/xuebo/Projects/Speciation/introgression/scanAe/raw_withbarleyAe/row_chr${chr}.withBarleyAe.vcf \
	--out /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr${chr}.withBarleyAe.vcf &
done
nohup sh getwithBarleyAe2.sh > nohup2.out 2>& 1 &

##因为至少加上两个outgroup,D lineage使用uratu和barley作为外类群，uratu要做scan，uratu往D基因组上做的mapping
/data1/home/yaozhou/data/ref/wheat/genome/outgroup/bwa/Tu 这是周姚的数据存放的位置
/data2/xuebo/Projects/Speciation/introgression/outgroup/Tu  这是204上面的位置
#################################################scan
for i in {1..42}
do
	grep -v "#" /data2/xuebo/Projects/Speciation/E3/chr${i}.all.vcf | cut -f 1-2,4-5 > chr${i}.allele.txt &
done
for i in {1..42}
do
	sed -i '1i Chr\tPos\tRef\tAlt' chr${i}.allele.txt &
done
for i in {1..42}
do
	grep -v "#" /data2/xuebo/Projects/Speciation/E3/chr${i}.all.vcf | cut -f 1-2 > chr${i}.pos.txt &
done
***** runningTu.sh 这个是D的外类群，所以只做D
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
    java -jar -Xms10g -Xmx20g /data2/xuebo/Projects/Speciation/hapScan/hapScanner2_jar/PlantGenetics.jar \
    /data2/xuebo/Projects/Speciation/introgression/scanTu/parameters_file/parameters_hapScannerTu_chr${chr}.txt > logchr${chr}.txt 
done
nohup sh runningTu.sh > nohupTu 2>& 1 &
***** 压缩
#!/bin/bash
for chr in {"005","006","011","012","017","018","023","024","029","030","035","036","041","042"}
do
	  i=$(echo $chr | sed 's/^0*//g')
	  #echo $i
	  bgzip -c /data2/xuebo/Projects/Speciation/introgression/scanTu/outchr${i}/VCF/chr${chr}.vcf > /data2/xuebo/Projects/Speciation/introgression/scanTu/scan_Tu_vcf/chr${i}.Tu.vcf.gz &
done
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
	tabix chr${chr}.Tu.vcf.gz &
done
****合并，查看行数  行数正确
grep -v "#" /data2/xuebo/Projects/Speciation/introgression/scanTu/outchr5/VCF/chr005.vcf | wc -l &
grep -v "#" /data2/xuebo/Projects/Speciation/introgression/scanTu/outchr11/VCF/chr011.vcf | wc -l &
*merge2.sh
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
		bcftools merge  -Oz /data2/xuebo/Projects/Speciation/tree/withBarley/AD_noblank/row_chr${chr}.noblank_withBarley.vcf.gz /data2/xuebo/Projects/Speciation/introgression/scanTu/scan_Tu_vcf/chr${chr}.Tu.vcf.gz \
    -o /data2/xuebo/Projects/Speciation/introgression/scanTu/raw_withbarleyTu/row_chr${chr}.withBarleyTu.vcf.gz &
done
nohup sh merge2.sh > nohupmerge2 2>& 1 &
zcat row_chr5.withBarleyTu.vcf.gz | grep -v "#" | wc -l
zcat row_chr6.withBarleyTu.vcf.gz | grep -v "#" | wc -l
*tbi.sh
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
  tabix row_chr${chr}.withBarleyTu.vcf.gz &
done
###把VCF文件里面Ae和barley是./.的去掉，保留Ae和barley都是有数值的
getwithBarleyTu.sh 
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/C26_getwithBarley_Ae_VCF.jar --file1 /data2/xuebo/Projects/Speciation/introgression/scanTu/raw_withbarleyTu/row_chr${chr}.withBarleyTu.vcf.gz \
	--out /data2/xuebo/Projects/Speciation/introgression/scanTu/withbarleyTu/chr${chr}.withBarleyTu.vcf &
	monitor java 20 60s
done
nohup sh getwithBarleyTu.sh > nohup.out 2>& 1 &

############################这是AB的加上的是Barley和Ae作为外类群 /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr${chr}.withBarleyAe.vcf
############################这是D的加上的是Barley和Tu作为外类群 /data2/xuebo/Projects/Speciation/introgression/scanTu/withbarleyTu/chr${chr}.withBarleyTu.vcf

############开始计算introgression  /data2/xuebo/Projects/Speciation/introgression
###先得到lineage的文件
**A
cat chr2.withBarleyAe.vcf chr7.withBarleyAe.vcf chr8.withBarleyAe.vcf chr13.withBarleyAe.vcf chr14.withBarleyAe.vcf chr19.withBarleyAe.vcf chr20.withBarleyAe.vcf chr25.withBarleyAe.vcf chr26.withBarleyAe.vcf chr31.withBarleyAe.vcf chr32.withBarleyAe.vcf chr37.withBarleyAe.vcf chr38.withBarleyAe.vcf | grep -v "#" > /data2/xuebo/Projects/Speciation/introgression/genoFile/Alineage/noA1.vcf
cat /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr1.withBarleyAe.vcf noA1.vcf > Alineage_withBarleyAe.vcf
bgzip -c Alineage_withBarleyAe.vcf > Alineage_withBarleyAe.vcf.gz &
tabix  Alineage_withBarleyAe.vcf.gz &
**B
cat chr4.withBarleyAe.vcf chr9.withBarleyAe.vcf chr10.withBarleyAe.vcf chr15.withBarleyAe.vcf chr16.withBarleyAe.vcf chr21.withBarleyAe.vcf chr22.withBarleyAe.vcf chr27.withBarleyAe.vcf chr28.withBarleyAe.vcf chr33.withBarleyAe.vcf chr34.withBarleyAe.vcf chr39.withBarleyAe.vcf chr40.withBarleyAe.vcf | grep -v "#" > /data2/xuebo/Projects/Speciation/introgression/genoFile/Blineage/noB1.vcf
cat /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr3.withBarleyAe.vcf noB1.vcf > Blineage_withBarleyAe.vcf
bgzip -c Blineage_withBarleyAe.vcf > Blineage_withBarleyAe.vcf.gz &
tabix  Blineage_withBarleyAe.vcf.gz &
**AB 
bcftools view -Oz -S /data2/xuebo/Projects/Speciation/introgression/groupforIntro/ABlineage_forfd.txt Alineage_withBarleyAe.vcf.gz --threads 10 -o ../ABlineage/lineageAB_A.vcf.gz &
bcftools view -Oz -S /data2/xuebo/Projects/Speciation/introgression/groupforIntro/ABlineage_forfd.txt Blineage_withBarleyAe.vcf.gz --threads 10 -o ../ABlineage/lineageAB_B.vcf.gz &
grep -v "#" lineageAB_B.vcf > lineageAB_B_nohead.vcf
cat lineageAB_A.vcf lineageAB_B_nohead.vcf > ABlineage_withBarley.vcf
bgzip -c ABlineage_withBarley.vcf > ABlineage_withBarley.vcf.gz &
tabix  ABlineage_withBarley.vcf.gz &
**D
cat chr6.withBarleyTu.vcf chr11.withBarleyTu.vcf chr12.withBarleyTu.vcf chr17.withBarleyTu.vcf chr18.withBarleyTu.vcf chr23.withBarleyTu.vcf chr24.withBarleyTu.vcf chr29.withBarleyTu.vcf chr30.withBarleyTu.vcf chr35.withBarleyTu.vcf chr36.withBarleyTu.vcf chr41.withBarleyTu.vcf chr42.withBarleyTu.vcf | grep -v "#" > /data2/xuebo/Projects/Speciation/introgression/genoFile/Dlineage/noD1.vcf
cat /data2/xuebo/Projects/Speciation/introgression/scanTu/withbarleyTu/chr5.withBarleyTu.vcf noD1.vcf > Dlineage_withBarleyTu.vcf
bgzip -c Dlineage_withBarleyTu.vcf > Dlineage_withBarleyTu.vcf.gz &
tabix  Dlineage_withBarleyTu.vcf.gz &

############################################数据转换，改格式 /data2/xuebo/Projects/Speciation/introgression/genoFile
WGS --model vcf --type toDstatistic --file /data2/xuebo/Projects/Speciation/introgression/genoFile/ABlineage/ABlineage_withBarley.vcf --out /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno &
bgzip -c AB.geno > AB.geno.gz
WGS --model vcf --type toDstatistic --file /data2/xuebo/Projects/Speciation/introgression/genoFile/Alineage/Alineage_withBarleyAe.vcf --out /data2/xuebo/Projects/Speciation/introgression/genoFile/A.geno &
bgzip -c A.geno > A.geno.gz &
WGS --model vcf --type toDstatistic --file /data2/xuebo/Projects/Speciation/introgression/genoFile/Blineage/Blineage_withBarleyAe.vcf --out /data2/xuebo/Projects/Speciation/introgression/genoFile/B.geno &
bgzip -c B.geno > B.geno.gz
WGS --model vcf --type toDstatistic --file /data2/xuebo/Projects/Speciation/introgression/genoFile/Dlineage/Dlineage_withBarleyTu.vcf --out /data2/xuebo/Projects/Speciation/introgression/genoFile/D.geno &
bgzip -c D.geno > D.geno.gz

############################################开始跑数据了
************************************************************************************************************run_emmer
#!/bin/bash
python /data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/Dstatistic/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g ../genoFile/AB.geno.gz -P1 North_wildemmer -P2 Domesticated_emmer  -P3 South_wildemmer -O Outgroup  \
-o AB_emmer1.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile ../groupforIntro/emmer1.txt --writeFailedWindows &
##在周姚的服务器上面做
python /data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/Dstatistic/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g ../genoFile/AB.geno.gz -P1 North_wildemmer -P2 Domesticated_emmer  -P3 South_wildemmer -O outgroup  \
-o AB_emmer1.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile ../groupforIntro/emmer1.txt --writeFailedWindows &
#
python /data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/Dstatistic/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g ../genoFile/AB.geno.gz -P1 Indian_dwarf_wheat -P2 Spelt  -P3 Domesticated_emmer -O outgroup  \
-o domemmer_spelt.csv.gz -w 100 --overlap 50 -T 80 -m 3 --popsFile ../groupforIntro/allABlineage.txt --writeFailedWindows &
python /data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/Dstatistic/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g ../genoFile/AB.geno.gz -P1 Indian_dwarf_wheat -P2 Xinjiang_wheat  -P3 Rivet_wheat -O outgroup  \
-o free_xinjiang.csv.gz -w 100 --overlap 50 -T 80 -m 3 --popsFile ../groupforIntro/allABlineage.txt --writeFailedWindows &
#
python /data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/Dstatistic/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g ../genoFile/AB.geno.gz -P1 EU_landrace -P2 Persian_wheat  -P3 EA_landrace -O outgroup  \
-o EAland_Persian.csv.gz -w 100 --overlap 50 -T 100 -m 3 --popsFile ../groupforIntro/allABlineage_landrace.txt --writeFailedWindows &

************************************************************************************************************Uratu_tatraploids
####之前跑的之一有问题，是python的版本问题
conda create -n py27 python=2.7
# /data2/xuebo/Projects/Speciation/introgression/fd/Uratu_tatraploids
#!/bin/bash
cat Info_Uratu_tatraploids_Pgroup.txt | while read line
do
	P1=$(echo $line | cut -f 1 -d$' ')
	P2=$(echo $line | cut -f 2 -d$' ')
	P3=$(echo $line | cut -f 3 -d$' ')
	P4=$(echo $line | cut -f 4 -d$' ')
	echo $P1
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/A.geno.gz \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${P4} -o fd_Uratu_${P2}.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/Info_Uratu_tatraploids.txt \
	--writeFailedWindows &
done
#把合适的fd挑选出来
ls -l |awk '{print $NF}' > list.txt
#!/usr/bin/Rscript
`getfd` <- function(file,name){
  f = read.csv(file,header = T)
  index = f$D>0 & f$fd < 1 & f$fd>0
  fdv = f[index,]$fd
  fdv = na.omit(fdv)
  fd = data.frame(ID=name,fd=fdv)
  return(fd)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Uratu_tatraploids")
fdvalue = NULL
file.index = c(1:9)
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1], as.character)
allfd = NULL
for(i in 1:length(file.index)){
  print(paste("Reading file ",info[file.index[i],1],sep=""))
  name = strsplit(info[file.index[i],1],split="\\.")
  name = unlist(name)[1]
  fd = getfd(info[file.index[i],1],name)
  allfd = rbind(allfd,fd)
}
write.table(allfd,"./usablefd/general_Uratu_tatraploids.txt",col.names = T,row.names = F,quote=F,sep="\t")
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Uratu_tatraploids")
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1],as.character)
ratio = data.frame(name = info[,1])
ratio$ratio = NA
for(i in 1:nrow(info)){
  print(paste("Current file is: ",info[i,1],sep=""))
  ratio[i,2] = fdRatio(info[i,1],T)
  # ratio_p[i,2] = fdRatioZtest(info[i,1])
}
write.table(ratio,"./usablefd/proportion_Uratu_tatraploids.txt",col.names=T,row.names=F,quote=F,sep="\t")

************************************************************************************************************Speltoides_tatraploids
# /data2/xuebo/Projects/Speciation/introgression/fd/Uratu_tatraploids
#!/bin/bash
cat Info_Speltoides_tatraploids_Pgroup.txt | while read line
do
	P1=$(echo $line | cut -f 1 -d$' ')
	P2=$(echo $line | cut -f 2 -d$' ')
	P3=$(echo $line | cut -f 3 -d$' ')
	P4=$(echo $line | cut -f 4 -d$' ')
	echo $P1
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/B.geno.gz \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${P4} -o fd_Speltoides_${P2}.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/Info_Speltoides_tatraploids.txt \
	--writeFailedWindows &
done
#把合适的fd挑选出来
ls -l |awk '{print $NF}' > list.txt
#!/usr/bin/Rscript
`getfd` <- function(file,name){
  f = read.csv(file,header = T)
  index = f$D>0 & f$fd < 1 & f$fd>0
  fdv = f[index,]$fd
  fdv = na.omit(fdv)
  fd = data.frame(ID=name,fd=fdv)
  return(fd)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Speltoides_tatraploids")
fdvalue = NULL
file.index = c(1:9)
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1], as.character)
allfd = NULL
for(i in 1:length(file.index)){
  print(paste("Reading file ",info[file.index[i],1],sep=""))
  name = strsplit(info[file.index[i],1],split="\\.")
  name = unlist(name)[1]
  fd = getfd(info[file.index[i],1],name)
  allfd = rbind(allfd,fd)
}
write.table(allfd,"./usablefd/general_Speltoides_tatraploids.txt",col.names = T,row.names = F,quote=F,sep="\t")
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Speltoides_tatraploids")
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1],as.character)
ratio = data.frame(name = info[,1])
ratio$ratio = NA
for(i in 1:nrow(info)){
  print(paste("Current file is: ",info[i,1],sep=""))
  ratio[i,2] = fdRatio(info[i,1],T)
  # ratio_p[i,2] = fdRatioZtest(info[i,1])
}
write.table(ratio,"./usablefd/proportion_Speltoides_tatraploids.txt",col.names=T,row.names=F,quote=F,sep="\t")

************************************************************************************************************Speltoides_domemmer
# /data2/xuebo/Projects/Speciation/introgression/fd/Speltoides_domemmer
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/B.geno.gz \
	-P1 Wild_emmer -P2 Domesticated_emmer -P3 Speltoides -O outgroup -o fd_Speltoides_domemmer.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/Info_Speltoides_domemmer.txt \
	--writeFailedWindows &
#把合适的fd挑选出来这一步没有写R，因为直接在R里面把全基因的图画一下即可
#算比例
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
ratio = data.frame(name = "Speltoides_domemmer")
ratio$ratio = NA
ratio[1,2] = fdRatio("fd_Speltoides_domemmer.csv.gz",T)
write.table(ratio,"./usablefd/proportion_Speltoides_domemmer.txt",col.names=T,row.names=F,quote=F,sep="\t")
##305843558-308909914的这个区域的introgression的位置是introgression比较高的位置--3M   300562860-327838165--20M
bcftools view -r 3:305843558-308909914 /data2/xuebo/Projects/Speciation/E3/beagle/chr3.beagle.vcf.gz -o fd_1B_3M.vcf
bcftools view -r 3:305843558-306578218 /data2/xuebo/Projects/Speciation/E3/beagle/chr3.beagle.vcf.gz -o fd_1B_1M.vcf

************************************************************************************************************Tatraploids_Hexa,算了indian dwarf 和yunan的P1，这里只展示indian的代码
# /data2/xuebo/Projects/Speciation/introgression/fd/Tatraploids_Hexa
#!/bin/bash
cat Info_Tatraploids_Hexa_Pgroup.txt | while read line
do
	P1=$(echo $line | cut -f 1 -d$' ')
	P2=$(echo $line | cut -f 2 -d$' ')
	P3=$(echo $line | cut -f 3 -d$' ')
	P4=$(echo $line | cut -f 4 -d$' ')
	echo $P1
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${P4} -o fd_Tatraploids_${P2}.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/Info_Tatraploids_Hexa.txt \
	--writeFailedWindows &
done
#把合适的fd挑选出来
ls -l |awk '{print $NF}' > list.txt
#!/usr/bin/Rscript
`getfd` <- function(file,name){
  f = read.csv(file,header = T)
  index = f$D>0 & f$fd < 1 & f$fd>0
  fdv = f[index,]$fd
  fdv = na.omit(fdv)
  fd = data.frame(ID=name,fd=fdv)
  return(fd)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Tatraploids_Hexa")
fdvalue = NULL
file.index = c(1:9)
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1], as.character)
allfd = NULL
for(i in 1:length(file.index)){
  print(paste("Reading file ",info[file.index[i],1],sep=""))
  name = strsplit(info[file.index[i],1],split="\\.")
  name = unlist(name)[1]
  fd = getfd(info[file.index[i],1],name)
  allfd = rbind(allfd,fd)
}
write.table(allfd,"./usablefd/general_Tatraploids_Hexa.txt",col.names = T,row.names = F,quote=F,sep="\t")
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Tatraploids_Hexa")
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1],as.character)
ratio = data.frame(name = info[,1])
ratio$ratio = NA
for(i in 1:nrow(info)){
  print(paste("Current file is: ",info[i,1],sep=""))
  ratio[i,2] = fdRatio(info[i,1],T)
  # ratio_p[i,2] = fdRatioZtest(info[i,1])
}
write.table(ratio,"./usablefd/proportion_Tatraploids_Hexa.txt",col.names=T,row.names=F,quote=F,sep="\t")

************************************************************************************************************Strangulata_Hexa
# /data2/xuebo/Projects/Speciation/introgression/fd/Strangulata_Hexa
#!/bin/bash
cat Info_Strangulata_Hexa_Pgroup.txt | while read line
do
	P1=$(echo $line | cut -f 1 -d$' ')
	P2=$(echo $line | cut -f 2 -d$' ')
	P3=$(echo $line | cut -f 3 -d$' ')
	P4=$(echo $line | cut -f 4 -d$' ')
	echo $P1
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/D.geno.gz \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${P4} -o fd_Strangulata_${P2}.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/Info_Strangulata_Hexa.txt \
	--writeFailedWindows &
done
#把合适的fd挑选出来
ls -l |awk '{print $NF}' > list.txt
#!/usr/bin/Rscript
`getfd` <- function(file,name){
  f = read.csv(file,header = T)
  index = f$D>0 & f$fd < 1 & f$fd>0
  fdv = f[index,]$fd
  fdv = na.omit(fdv)
  fd = data.frame(ID=name,fd=fdv)
  return(fd)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Strangulata_Hexa")
fdvalue = NULL
file.index = c(1:9)
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1], as.character)
allfd = NULL
for(i in 1:length(file.index)){
  print(paste("Reading file ",info[file.index[i],1],sep=""))
  name = strsplit(info[file.index[i],1],split="\\.")
  name = unlist(name)[1]
  fd = getfd(info[file.index[i],1],name)
  allfd = rbind(allfd,fd)
}
write.table(allfd,"./usablefd/general_Strangulata_Hexa.txt",col.names = T,row.names = F,quote=F,sep="\t")
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Strangulata_Hexa")
info = read.table("list.txt",head=F)
info[,1] = sapply(info[,1],as.character)
ratio = data.frame(name = info[,1])
ratio$ratio = NA
for(i in 1:nrow(info)){
  print(paste("Current file is: ",info[i,1],sep=""))
  ratio[i,2] = fdRatio(info[i,1],T)
  # ratio_p[i,2] = fdRatioZtest(info[i,1])
}
write.table(ratio,"./usablefd/proportion_Strangulata_Hexa.txt",col.names=T,row.names=F,quote=F,sep="\t")


************************************************************************************************************Spelt
# /data2/xuebo/Projects/Speciation/introgression/fd/Spelt
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Indian_dwarf_wheat -P2 Spelt -P3 Domesticated_emmer -O outgroup -o fd_domemmer_spelt.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Cultivar -P2 Spelt -P3 Landrace -O outgroup -o fd_landrace_spelt.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/spelt")
ratio = data.frame(name = "spelt")
ratio$name = NA
ratio$ratio = NA
ratio[1,1] = "domemmer_spelt"
ratio[2,1] = "landrace_spelt"
ratio[1,2] = fdRatio("fd_domemmer_spelt.csv.gz",T)
ratio[2,2] = fdRatio("fd_landrace_spelt.csv.gz",T)
write.table(ratio,"./ratio_spelt.txt",col.names=T,row.names=F,quote=F,sep="\t")


************************************************************************************************************Spelt
# /data2/xuebo/Projects/Speciation/introgression/fd/Spelt
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Indian_dwarf_wheat -P2 Spelt -P3 Domesticated_emmer -O outgroup -o fd_domemmer_spelt.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Cultivar -P2 Spelt -P3 Landrace -O outgroup -o fd_landrace_spelt.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/spelt")
ratio = data.frame(name = "spelt")
ratio$name = NA
ratio$ratio = NA
ratio[1,1] = "domemmer_spelt"
ratio[2,1] = "landrace_spelt"
ratio[1,2] = fdRatio("fd_domemmer_spelt.csv.gz",T)
ratio[2,2] = fdRatio("fd_landrace_spelt.csv.gz",T)
write.table(ratio,"./ratio_spelt.txt",col.names=T,row.names=F,quote=F,sep="\t")


************************************************************************************************************Macha
# /data2/xuebo/Projects/Speciation/introgression/fd/Macha
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Indian_dwarf_wheat -P2 Macha -P3 Domesticated_emmer -O outgroup -o fd_domemmer_Macha.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Cultivar -P2 Macha -P3 Landrace -O outgroup -o fd_landrace_Macha.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Macha")
ratio = data.frame(name = "Macha")
ratio$name = NA
ratio$ratio = NA
ratio[1,1] = "domemmer_Macha"
ratio[2,1] = "landrace_Macha"
ratio[1,2] = fdRatio("fd_domemmer_Macha.csv.gz",T)
ratio[2,2] = fdRatio("fd_landrace_Macha.csv.gz",T)
write.table(ratio,"./ratio_Macha.txt",col.names=T,row.names=F,quote=F,sep="\t")


************************************************************************************************************Xinjiang_wheat
# /data2/xuebo/Projects/Speciation/introgression/fd/Xinjiang_wheat
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Indian_dwarf_wheat -P2 Xinjiang_wheat -P3 Polish_wheat -O outgroup -o fd_Polish_wheat_Xinjiang_wheat.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Cultivar -P2 Xinjiang_wheat -P3 Landrace -O outgroup -o fd_landrace_Xinjiang_wheat.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Xinjiang_wheat")
ratio = data.frame(name = "Xinjiang_wheat")
ratio$name = NA
ratio$ratio = NA
ratio[1,1] = "Polish_wheat_Xinjiang_wheat"
ratio[2,1] = "landrace_Xinjiang_wheat"
ratio[1,2] = fdRatio("fd_Polish_wheat_Xinjiang_wheat.csv.gz",T)
ratio[2,2] = fdRatio("fd_landrace_Xinjiang_wheat.csv.gz",T)
write.table(ratio,"./ratio_Xinjiang_wheat.txt",col.names=T,row.names=F,quote=F,sep="\t")

************************************************************************************************************Persian_wheat
# /data2/xuebo/Projects/Speciation/introgression/fd/Persian_wheat
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Indian_dwarf_wheat -P2 Persian_wheat -P3 Rivet_wheat -O outgroup -o fd_Rivet_wheat_Persian_wheat.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#!/bin/bash
	python2 /data2/xuebo/Projects/Speciation/introgression/genomics_general/ABBABABAwindows.py  --windType sites -f phased -g /data2/xuebo/Projects/Speciation/introgression/genoFile/AB.geno.gz \
	-P1 Cultivar -P2 Persian_wheat -P3 Landrace -O outgroup -o fd_landrace_Persian_wheat.csv.gz -w 100 --overlap 50 -T 30 -m 3 --popsFile /data2/xuebo/Projects/Speciation/introgression/groupforIntro/allABlineage.txt \
	--writeFailedWindows &
#算比例
#!/usr/bin/Rscript
`fdRatio` <- function(file,csv=F){
  if(csv){
    res = read.csv(file,header = T)
  }else{
    res = read.table(file,head=T)
  }
  chr = unique(res$scaffold)
  sum_fd = 0
  sum_all = 0
  for (i in chr){
    index = res$scaffold == i
    dat = res[index,]
    dat2 = dat
    dat2[1:(nrow(dat)-1),3] = dat[2:nrow(dat),2]
    index = dat2$D > 0 & dat2$fd < 1 & dat2$fd > 0
    sum = (dat2[index,]$end - dat2[index,]$start)*dat2[index,]$fd
    sum_fd = sum(sum,na.rm=T) + sum_fd
    sum_all = sum((dat2$end - dat2$start),na.rm=T) + sum_all
  }
  return(sum_fd/sum_all)
}
setwd("/data2/xuebo/Projects/Speciation/introgression/fd/Persian_wheat")
ratio = data.frame(name = "Persian_wheat")
ratio$name = NA
ratio$ratio = NA
ratio[1,1] = "Rivet_wheat_Persian_wheat"
ratio[2,1] = "landrace_Persian_wheat"
ratio[1,2] = fdRatio("fd_Rivet_wheat_Persian_wheat.csv.gz",T)
ratio[2,2] = fdRatio("fd_landrace_Persian_wheat.csv.gz",T)
write.table(ratio,"./ratio_Persian_wheat.txt",col.names=T,row.names=F,quote=F,sep="\t")


































