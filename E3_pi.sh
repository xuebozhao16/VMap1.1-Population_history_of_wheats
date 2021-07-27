############ Blineage的pi情况
//WGS --model diversity --type bedPi --size 10000 --file lineageA.1M.group_cultivar_small.sites.pi --file2 /data1/home/yaozhou/Projects/EVO/data/merge/lineage/bed/all.bed2 --out ./Pi10k/lineageA.10k.cultivar_small.bed.pi &
***site pi
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
	vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz  --keep /data2/xuebo/Projects/Speciation/group/SS_taxa.txt --site-pi  --out chr${chr}.SS  &
done
***计算bed里面有多少SNP  /data2/xuebo/Projects/Speciation/E3/
#!/bin/bash
for chr in {1..42}
do
  WGS --model depth --type toBed --file /data2/xuebo/Projects/Speciation/E3/synteny_site/chr${chr}.txt --out /data2/xuebo/Projects/Speciation/E3/SNPbed_10K/chr${chr}_10000.bed --size 10000 &
done
#!/bin/bash
for chr in {1..42}
do
  WGS --model depth --type toBed --file /data2/xuebo/Projects/Speciation/E3/synteny_site/chr${chr}.txt --out /data2/xuebo/Projects/Speciation/E3/SNPbed_100K/chr${chr}_100K.bed --size 100000 &
done
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
  WGS --model depth --type toBed --file /data2/xuebo/Projects/Speciation/E3/synteny_all_B/chr${chr}.txt --out /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed --size 1000000 &
done
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
  WGS --model depth --type toBed --file /data2/xuebo/Projects/Speciation/E3/synteny_siteV11/chr${chr}.txt --out /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed --size 1000000 &
done
***计算bed pi 
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do 
  WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/SS/chr${chr}.SS.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
  --out chr${chr}.SS.bed.pi &
done
cat chr* | sort -k1,1 -k2,2n > B.group_speltoides.bed.pi


######################计算wild emmer的南北分群的多样性的差异
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
	vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz  --keep /data2/xuebo/Projects/Speciation/group/wildemmer/South_wildemmer.txt --site-pi  --out chr${chr}.South_wildemmer  &
	vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz  --keep /data2/xuebo/Projects/Speciation/group/wildemmer/North_wildemmer.txt --site-pi  --out chr${chr}.North_wildemmer  &
	vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz  --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Domesticated_emmer.txt --site-pi  --out chr${chr}.Domesticated_emmer  &
done
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do 
  WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.South_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
  --out chr${chr}.South_wildemmer.bed.pi &
  WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.North_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
  --out chr${chr}.North_wildemmer.bed.pi &
  WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.Domesticated_emmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
  --out chr${chr}.Domesticated_emmer.bed.pi &
done
cat chr*.South_wildemmer.bed.pi | sort -k1,1 -k2,2n > South_wildemmer.bed.pi
cat chr*.North_wildemmer.bed.pi | sort -k1,1 -k2,2n > North_wildemmer.bed.pi
cat chr*.Domesticated_emmer.bed.pi | sort -k1,1 -k2,2n > Domesticated_emmer.bed.pi
*****分别计算南北方的emmer的多样性
#!/bin/bash
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do 
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.South_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.South_wildemmer_A.bed.pi &
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.North_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.North_wildemmer.bed_A.pi &
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.Domesticated_emmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.Domesticated_emmer_A.bed.pi &
done
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do 
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.South_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.South_wildemmer_B.bed.pi &
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.North_wildemmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.North_wildemmer.bed_B.pi &
	WGS --chr ${chr} --model diversity --type bedPi --file /data2/xuebo/Projects/Speciation/diversity/emmer/chr${chr}.Domesticated_emmer.sites.pi --file2 /data2/xuebo/Projects/Speciation/E3/SNPbed_1M/chr${chr}_1M.bed \
	--out chr${chr}.Domesticated_emmer_B.bed.pi &
done
cat chr*.South_wildemmer_A.bed.pi | sort -k1,1 -k2,2n > South_wildemmer_A.bed.pi
cat chr*.South_wildemmer_B.bed.pi | sort -k1,1 -k2,2n > South_wildemmer_B.bed.pi
cat chr*.North_wildemmer.bed_A.pi | sort -k1,1 -k2,2n > North_wildemmer_A.bed.pi
cat chr*.North_wildemmer.bed_B.pi | sort -k1,1 -k2,2n > North_wildemmer_B.bed.pi
cat chr*.Domesticated_emmer_A.bed.pi | sort -k1,1 -k2,2n > Domesticated_emmer_A.bed.pi
cat chr*.Domesticated_emmer_B.bed.pi | sort -k1,1 -k2,2n > Domesticated_emmer_B.bed.pi
#############PCA分析 
conda install -c bioconda plink //这个不行
wget -c http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip
unzip plink_linux_x86_64_20190617.zip
PATH=$PATH:/data1/home/xuebo/software/
source ~/.bashrc  //但是这样是把整个software加进环境变量，这样不好
可以把plink这个可执行程序加进已经有的环境环境变量里面   PATH=$PATH:/data1/home/xuebo/software/bedtools2/bin
mv plink /data1/home/xuebo/software/bedtools2/bin
source ~/.bashrc //这样就好啦
先用画树的with barley的100K的SNP做一个测试版本
***首先把barley去掉
vcftools --gzvcf /data2/xuebo/Projects/Speciation/tree/random/Alineage/Atree.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/Alineage_taxa.txt --recode --stdout > /data2/xuebo/Projects/Speciation/pca/random100K/Atree100k.vcf &
vcftools --gzvcf /data2/xuebo/Projects/Speciation/tree/random/Blineage/Btree.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/Blineage_taxa.txt --recode --stdout > /data2/xuebo/Projects/Speciation/pca/random100K/Btree100k.vcf &
vcftools --gzvcf /data2/xuebo/Projects/Speciation/tree/random/Dlineage/Dtree.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/Dlineage_taxa.txt --recode --stdout > /data2/xuebo/Projects/Speciation/pca/random100K/Dtree100k.vcf &

bgzip -c Atree100k.vcf > Atree100k.vcf.gz &
bgzip -c Btree100k.vcf > Btree100k.vcf.gz &
bgzip -c Dtree100k.vcf > Dtree100k.vcf.gz &

tabix Atree100k.vcf.gz
tabix Btree100k.vcf.gz
tabix Dtree100k.vcf.gz

vcftools --gzvcf ../Atree100k.vcf.gz --plink --out A.lineage &
vcftools --gzvcf ../Btree100k.vcf.gz --plink --out B.lineage &
vcftools --gzvcf ../Dtree100k.vcf.gz --plink --out D.lineage &

plink --file A.lineage --pca 100 header tabs --chr-set 42 --out A &
plink --file B.lineage --pca 100 header tabs --chr-set 42 --out B &
plink --file D.lineage --pca 100 header tabs --chr-set 42 --out D &

plink --file AexontreeIN  --pca 100 header tabs --chr-set 42 --out AAexontreeIN &

##现在算的是ABD里面的D里面的PCA /data2/xuebo/Projects/Speciation/pca/D_inABD
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
    vcftools --vcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf --keep /data2/xuebo/Projects/Speciation/group/AABBDD_taxa.txt \
    --maf 0.00001  --recode --stdout > /data2/xuebo/Projects/Speciation/pca/D_inABD/vcffile/chr${chr}.ABDtaxa.vcf &
done
#!/bin/bash
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
	WGS --model vcf --type overlap --file /data2/xuebo/Projects/Speciation/pca/D_inABD/vcffile/chr${chr}.ABDtaxa.vcf \
	--file2 /data2/xuebo/Projects/Speciation/gff/chr${chr}.exon --out /data2/xuebo/Projects/Speciation/pca/D_inABD/vcffile/chr${chr}.genictree.vcf &
done
vcf-concat chr5.genictree.vcf chr6.genictree.vcf chr11.genictree.vcf chr12.genictree.vcf chr17.genictree.vcf chr18.genictree.vcf chr23.genictree.vcf chr24.genictree.vcf chr29.genictree.vcf chr30.genictree.vcf chr35.genictree.vcf chr36.genictree.vcf chr41.genictree.vcf chr42.genictree.vcf > Dlieange_ABDtaxa.vcf
88,256
bgzip  Dlieange_ABDtaxa.vcf
tabix Dlieange_ABDtaxa.vcf.gz
vcftools --gzvcf Dlieange_ABDtaxa.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/AABBDD_taxa_NoSynthetic.txt --maf 0.00001  --recode --stdout > Dlieange_ABDtaxa2.vcf
78,164
bgzip  Dlieange_ABDtaxa2.vcf
tabix Dlieange_ABDtaxa2.vcf.gz
vcftools --gzvcf Dlieange_ABDtaxa2.vcf.gz --plink --out Dlieange_ABDtaxa &
plink --file Dlieange_ABDtaxa --pca 100 header tabs --chr-set 42 --out D.ABDtaxa &

###########################现在是使用画树的那套SNP做PCA,计算MDS 
#/data2/xuebo/Projects/Speciation/pca/E3_tree_vcf
cp AexontreeIN.vcf.gz* /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome
cp BexontreeIN.vcf.gz* /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome
vcftools --gzvcf AexontreeIN.vcf.gz --plink --out AexontreeIN 
plink --file AexontreeIN --distance-matrix --chr-set 42 --out A_IN_matrix
plink --file AexontreeIN --distance 1-ibs --chr-set 42 --out A_IN_matrix ##这个是角阵
plink --file AexontreeIN --distance square '1-ibs' --chr-set 42 --out A_IN_matrix #方阵
vcftools --gzvcf BexontreeIN.vcf.gz --plink --out BexontreeIN 
plink --file BexontreeIN --distance square '1-ibs' --chr-set 42 --out B_IN_matrix #方阵
**这个B的结果和A的太像了，现在做一个check /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome/checkB
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
    vcftools --vcf /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf --keep /data2/xuebo/Projects/Speciation/group/Alineage_taxa_inner.txt \
    --maf 0.00001  --recode --stdout | bgzip  > /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome/checkB/chr${chr}.checkBtaxa.vcf.gz &
done
#!/bin/bash
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
	WGS --model vcf --type overlap --file /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome/checkB/vcf/chr${chr}.checkBtaxa.vcf.gz \
	--file2 /data2/xuebo/Projects/Speciation/gff/chr${chr}.exon --out /data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/ABgenome/checkB/Exon_chr${chr}.vcf &
done
vcf-concat Exon_chr3.vcf Exon_chr4.vcf Exon_chr9.vcf Exon_chr10.vcf Exon_chr15.vcf Exon_chr16.vcf Exon_chr21.vcf Exon_chr22.vcf Exon_chr27.vcf Exon_chr28.vcf Exon_chr33.vcf Exon_chr34.vcf Exon_chr39.vcf Exon_chr40.vcf > Exon_InB.vcf
212,386
vcftools --gzvcf Exon_InB.vcf.gz --plink --out BexontreeIN 
plink --file BexontreeIN --distance square '1-ibs' --chr-set 42 --out B_IN_matrix
#/data2/xuebo/Projects/Speciation/pca/E3_tree_vcf/D/IN_D
cp /data2/xuebo/Projects/Speciation/pca/D_inABD/vcffile/Dlieange_ABDtaxa2.vcf.gz* ./
vcftools --gzvcf Dlieange_ABDtaxa2.vcf.gz --plink --out DexontreeIN 
plink --file DexontreeIN --distance square '1-ibs' --chr-set 42 --out D_IN_matrix


###########################################################################sNMF
*mkdir.sh
#!/bin/bash
for i in {1..100}
do
	mkdir rep$i
done
*vcf2geno
vcf2geno /data2/xuebo/Projects/Speciation/tree/innertree/Alineage/Agenictree_inner.vcf   Ainner.geno
vcf2geno /data2/xuebo/Projects/Speciation/tree/innertree/Blineage/Bexontree_withBarleyNoS.vcf  Binner.geno
vcf2geno /data2/xuebo/Projects/Speciation/tree/innertree/Dlineage_maf/maf01/chr.Dsubgenome_300k_maf01.recode.vcf  Dinner.geno
*cp.sh
#!/bin/bash
for i in {1..100}
do
ln Ainner.geno ./rep$i/ &
ln Binner.geno ./rep$i/ &
ln Dinner.geno ./rep$i/ &
done
*run.sh
#!/bin/bash
for i in {2..10}
do
        count=$( ps -e | grep "sNMF"| wc -l )
        echo $count
        while [ $count -gt 0 ]
        do
                count=$( ps -e | grep "sNMF"| wc -l )
                sleep 1s
        done
        for j in {1..100}
        do
                cd rep$j
                sNMF -x Ainner.geno -K $i -c > A.$i.log &
                sNMF -x Binner.geno -K $i -c > B.$i.log &
                sNMF -x Dinner.geno -K $i -c > D.$i.log &
                cd ..
        done
done
*Q2CLUMPP.sh 
#!/bin/bash
for i in {1..100}
do
        for j in {2..10}
        do
                WGS --model structure --type Q2CLUMPP --file ./rep$i/Ainner.${j}.Q --out ./rep$i/A.${j}.Q.m
                WGS --model structure --type Q2CLUMPP --file ./rep$i/Binner.${j}.Q --out ./rep$i/B.${j}.Q.m
                WGS --model structure --type Q2CLUMPP --file ./rep$i/Dinner.${j}.Q --out ./rep$i/D.${j}.Q.m
        done
done
*****/data2/xuebo/Projects/Speciation/sNMF/inner/CLUMPP
*mkdir.sh
for i in {2..10}
do
	mkdir K$i
done
*cat.sh
#!/bin/bash
for k in {2..10}
do
        cat ../rep1/A.$k.Q.m > ./K$k/A.Q
        cat ../rep1/B.$k.Q.m > ./K$k/B.Q
        cat ../rep1/D.$k.Q.m > ./K$k/D.Q
        for i in {2..100}
        do
                cat ../rep$i/A.$k.Q.m >> ./K$k/A.Q
                cat ../rep$i/B.$k.Q.m >> ./K$k/B.Q
                cat ../rep$i/D.$k.Q.m >> ./K$k/D.Q
        done
done
*getPara.sh
#!/bin/bash
for i in {3..10}
do
	java -jar ~/software/WGS.jar --model GenerateScripts --type model19 --file ./K2/paramfileA --out ./K$i/paramfileA --suffix $i
	java -jar ~/software/WGS.jar --model GenerateScripts --type model19 --file ./K2/paramfileB --out ./K$i/paramfileB --suffix $i
	java -jar ~/software/WGS.jar --model GenerateScripts --type model19 --file ./K2/paramfileD --out ./K$i/paramfileD --suffix $i
done
run.sh
#!/bin/bash
for i in {2..10}
do
        cd K$i
        CLUMPP paramfileA > A.log &
        CLUMPP paramfileB > B.log &
        CLUMPP paramfileD > D.log &
        cd ..
done
*CLUMPP2R.sh
#!/bin/bash
for i in {2..10}
do
WGS --model structure --type CLUMPP2R --file ./K$i/A.outfile --out ./result/A.${i}.Q &
WGS --model structure --type CLUMPP2R --file ./K$i/B.outfile --out ./result/B.${i}.Q &
WGS --model structure --type CLUMPP2R --file ./K$i/D.outfile --out ./result/D.${i}.Q &
done
*addlabel.sh
#!/bin/bash
for i in {2..10}
do
        paste /data2/xuebo/Projects/Speciation/group/Alineage_taxa_inner.txt A.$i.Q > labelA.$i.Q
        paste /data2/xuebo/Projects/Speciation/group/Blineage_taxa_inner.txt B.$i.Q > labelB.$i.Q
        paste /data2/xuebo/Projects/Speciation/group/Dlineage_taxa_inner.txt D.$i.Q > labelD.$i.Q
done



##TajimasD计算
#!/bin/bash
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Alineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Landrace.txt --TajimaD 1000000  --out lineageA.1M.group_Landrace  & 
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Blineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Landrace.txt --TajimaD 1000000  --out lineageB.1M.group_Landrace  & 
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Dlineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Landrace.txt --TajimaD 1000000  --out lineageD.1M.group_Landrace  & 
#!/bin/bash
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Alineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Cultivar.txt --TajimaD 1000000  --out lineageA.1M.group_Cultivar  & 
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Blineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Cultivar.txt --TajimaD 1000000  --out lineageB.1M.group_Cultivar  & 
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Dlineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Cultivar.txt --TajimaD 1000000  --out lineageD.1M.group_Cultivar  & 
#!/bin/bash
vcftools --gzvcf /data2/xuebo/Projects/Speciation/E3/Dlineage.vcf.gz --keep /data2/xuebo/Projects/Speciation/group/subspecies/sub_Strangulata.txt --TajimaD 1000000  --out lineageD.1M.group_Strangulata  & 










