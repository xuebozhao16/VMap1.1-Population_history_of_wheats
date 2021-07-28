##软件位置 data1/home/xuebo/software/fsc26_linux64/fsc26
##工作环境 203 /data1/home/xuebo/Projects/Speciation/fsc
##说明书 https://speciationgenomics.github.io/fastsimcoal2/  &&  http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf
https://github.com/speciationgenomics/presentations/blob/master/Demographic_modeling_JM.pdf
https://github.com/isaacovercast/easySFS
https://github.com/speciationgenomics/scripts
http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/inferring-demography/

#####这里面用的是加上外类群的SNP 这是在204上面的位置
这是AB的加上的是Barley和Ae作为外类群 /data2/xuebo/Projects/Speciation/introgression/scanAe/withbarleyAe/chr${chr}.withBarleyAe.vcf
这是D的加上的是Barley和Tu作为外类群 /data2/xuebo/Projects/Speciation/introgression/scanTu/withbarleyTu/chr${chr}.withBarleyTu.vcf
/data2/xuebo/Projects/Speciation/introgression/genoFile/Alineage/Alineage_withBarleyAe.vcf.gz
/data2/xuebo/Projects/Speciation/introgression/genoFile/Blineage/Blineage_withBarleyAe.vcf.gz
/data2/xuebo/Projects/Speciation/introgression/genoFile/Dlineage/Dlineage_withBarleyTu.vcf.gz

## /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS 这是以上文件在203上面的路径

###首先生成2D的SFS文件 工作路径/data1/home/xuebo/Projects/Speciation/fsc/2D_SFS
##说明文档 https://github.com/isaacovercast/easySFS
/data1/home/xuebo/software/easySFS/easySFS.py
***先看一下能不能运行
/data1/home/xuebo/software/easyåSFS/easySFS.py -i wcs_1200.vcf -p wcs_pops.txt --preview -a
/data1/home/xuebo/software/easySFS/easySFS.py -i wcs_1200.vcf -p wcs_pops.txt -a --proj=74,186
***现在是做D lineage
/data1/home/xuebo/software/easySFS/easySFS.py -i /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/Dlineage_withBarleyTu.vcf.gz -p D_pops.txt  --preview -a
/data1/home/xuebo/software/easySFS/easySFS.py -i /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/Dlineage_withBarleyTu.vcf.gz -p D_pops.txt -a --proj=40,40 -o /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/Dlineage/outputSFS
##现在准备其他的两个文件
**********************************************************************************先看个例子 /data1/home/xuebo/Projects/Speciation/fsc/example/pop2_mig_example
#!/bin/bash
PREFIX="2PopDivMigr20Mb"
/data1/home/xuebo/software/fsc26_linux64/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -d -0 -C 1 -n 10000 -L 40 -s 0 -M -c80
***********************现在是尝试不同的循环
#!/bin/bash
PREFIX="2PopDivMigr20Mb"
for i in {1..10}
do
	mkdir run$i
	cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}_jointDAFpop1_0.obs run$i"/"
	cd run$i
	/data1/home/xuebo/software/fsc26_linux64/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -d -0 -C 10 -n 10000 -L 40 -s0 -M -q -c80
	cd ..
done
************************现在是找到合适的那个
#!/bin/bash
PREFIX="2PopDivMigr20Mb"
cat run{1..10}/${PREFIX}/${PREFIX}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$8}' | sort -k 2  ##现在发现10个是一样的
************************现在是找到合适的那个，并且复制到一个叫做bestrun的文件里面 fsc-selectbestrun.sh
sh fsc-selectbestrun.sh
********************有5种不同的基因流模型，Different gene flow matrices这个模型看起来更合适，但是还是要用AIC表示模型的运行可能性
cd bestrun/
/data1/home/xuebo/Projects/Speciation/fsc/Rcode/calculateAIC.r 2PopDivMigr20Mb
********************可视化
cd bestrun/
/data1/home/xuebo/Projects/Speciation/fsc/Rcode/SFStools.r -t print2D -i 2PopDivMigr20Mb
Rscript /data1/home/xuebo/Projects/Speciation/fsc/Rcode/plotModel.r -p 2PopDivMigr20Mb -l NPOP1,NPOP2
**********************************************************************************现在做D lineage /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/Dlineage/
************************现在是找到合适的那个 分别建立了5个文件夹
************************检测5种模型是不是合适

###############现在是做AB lineage的，产生AB的obs文件
#wild emmer 0
#dom emmer 1
#free-threshing 2
#landrace 3
/data1/home/xuebo/software/easySFS/easySFS.py -i /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage_withBarley.vcf.gz -p AB_pops.txt -a --proj=40,40,40,40 -o /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/outputSFS 
cp ABlineage_withBarley_jointMAFpop1_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop1_0.obs 
cp ABlineage_withBarley_jointMAFpop2_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop2_0.obs 
cp ABlineage_withBarley_jointMAFpop3_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop3_0.obs 
cp ABlineage_withBarley_jointMAFpop2_1.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop2_1.obs 
cp ABlineage_withBarley_jointMAFpop3_1.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop3_1.obs 
cp ABlineage_withBarley_jointMAFpop3_2.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_jointMAFpop3_2.obs 
cp ABlineage_withBarley_MSFS.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_no/ABlineage_PopDiv_no_MSFS.obs

#!/bin/bash
PREFIX="ongoing"
cp ABlineage_withBarley_jointMAFpop1_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop1_0.obs 
cp ABlineage_withBarley_jointMAFpop2_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop2_0.obs 
cp ABlineage_withBarley_jointMAFpop3_0.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop3_0.obs 
cp ABlineage_withBarley_jointMAFpop2_1.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop2_1.obs 
cp ABlineage_withBarley_jointMAFpop3_1.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop3_1.obs 
cp ABlineage_withBarley_jointMAFpop3_2.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_jointMAFpop3_2.obs 
cp ABlineage_withBarley_MSFS.obs /data1/home/xuebo/Projects/Speciation/fsc/2D_SFS/ABlineage/model_${PREFIX}/ABlineage_PopDiv_${PREFIX}_MSFS.obs






###########################################################这是AB的,用的是最合适的那个模型
#参考文档 https://github.com/h-e-g/evoceania/blob/main/Fastsimcoal2_inputs/1-Baseline/Baseline_EUR_PNG_ASN.tpl
#### ABlineage_PopDiv_no.tpl
//Parameters for the coalescence simulation program : fastsimcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
WE
DE
FT
LAN
//Haploid samples sizes 
80
80
80
80
//Growth rates: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
MIG1 MIG2 MIG3 0
MIG4 MIG5 0 0
MIG6 0 0 0
0 0 0 0
//Migration matrix 1
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3 historical event
11000 1 0 1 RSANC 0 0
9500 2 1 0 RSANC 0 0
8000 3 2 0 RSANC 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 6.5e-9 OUTEXP

#### ABlineage_PopDiv_no.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all Ns are in number of haploid individuals
1  WE		logunif	100	 100000	 output
1  DE		logunif	100	 100000	 output
1  FT		logunif	100	 100000	 output
1  LAN		logunif	100	 100000	 output
0  RSANC	logunif	0.1	 100	 output
0  MIG1      	logunif	0.000001 0.001   output
0  MIG2      	logunif	0.000001 0.001   output
0  MIG3         logunif 0.000001 0.001   output
0  MIG4         logunif 0.000001 0.001   output
0  MIG5         logunif 0.000001 0.001   output
0  MIG6         logunif 0.000001 0.001   output
0  TPROP	logunif 0.0001	 0.5	 output

[RULES]

[C	OMPLEX PARAMETERS]


#####AB lienage
/data1/home/xuebo/software/fsc26_linux64/fsc26 -t ABlineage_PopDiv_no.tpl -e ABlineage_PopDiv_no.est -m -0 -C 1 -n 10 -L 40 -s 0 -M -c80
###############running_rep.sh
#!/bin/bash
PREFIX="ABlineage_PopDiv_no"
for i in {1..10}
do
	mkdir run$i
	cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}_*.obs run$i"/"
	cd run$i
	/data1/home/xuebo/software/fsc26_linux64/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -0 -C 10 -n 10000 -L 40 -s0 -M -q -c120
	cd ..
done

######get bestrun
*****
#!/bin/bash
PREFIX="ABlineage_PopDiv_no"
cat run{1..10}/${PREFIX}/${PREFIX}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$5"\t"$6}' | sort -k 2
*****
sh /data1/home/xuebo/Projects/Speciation/fsc/Rcode/fsc-selectbestrun.sh
*****AIC
cd bestrun/
/data1/home/xuebo/Projects/Speciation/fsc/Rcode/calculateAIC.r ABlineage_PopDiv_no
*****
cd bestrun/
Rscript /data1/home/xuebo/Projects/Speciation/fsc/Rcode/plotModel_early.r -p ABlineage_PopDiv_no -l WE,DE,FT,LAN





























