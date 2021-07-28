#####SwiftOrtho https://github.com/Rinoahu/SwiftOrtho
##得到不同的基因组的pepwe文件
*****测试
seqkit sort -n Triticum_aestivum.IWGSC.pep.all.fa.gz > Triticum_aestivum.IWGSC.pep.allsort.fa
grep -n "TraesCS1B" Triticum_aestivum.IWGSC.pep.allsort.fa | head -n 1  46536
head -n 46535 Triticum_aestivum.IWGSC.pep.allsort.fa > TaA.pep.fa
grep -n "TraesCS1D" Triticum_aestivum.IWGSC.pep.allsort.fa | head -n 1  98489
sed -n '46536,98488p' Triticum_aestivum.IWGSC.pep.allsort.fa > TaB.pep.fa
python /data1/home/xuebo/software/SwiftOrtho/bin/find_hit.py -p blastp -i /data1/home/xuebo/wheat/othergenome/urartu/Triticum_urartu.ASM34745v1.pep.all.fa.gz  -d /data1/home/xuebo/wheat/Ref/Triticum_aestivum.IWGSC.pep.all.fa.gz -o input.fsa.sc -e 1e-5 -s 111111 -a 80
python /data1/home/xuebo/software/SwiftOrtho/bin/find_orth.py -i Triticum_urartu.ASM34745v1.pep.all.fa.gz -c 0.5 -y 0 > input.fsa.sc.orth
zcat Triticum_aestivum.IWGSC.pep.all.fa.gz | awk -F " " '{print$1}' | head -n 2
zcat Triticum_aestivum.IWGSC.pep.all.fa.gz | cut -d '.' -f 1 | head -n 2 这个也可以用
zcat Triticum_aestivum.IWGSC.pep.all.fa.gz | awk '/>/ {$0 =substr($1,1,length($1)-2)} 1' > Triticum_aestivum.IWGSC.pep.allname.fa
cat Triticum_aestivum.IWGSC.pep.allname.fa|grep 'TraesCS4B' | awk '/>/ {if(substr($1,10,1)=="B"){$0=">TaB|"substr($1,2,length($1))}} 1' | head -n 2
*cat Triticum_aestivum.IWGSC.pep.allname.fa | awk '/>/ {if(substr($1,10,1)=="A"){$0=">TaA|"substr($1,2,length($1))}} 1' > Triticum_aestivum.IWGSC.pep.allnametmpA.fa
cat Triticum_aestivum.IWGSC.pep.allname.fa | awk '/>/ {if(substr($1,10,1)=="A"){$0=">TaA|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,10,1)=="B"){$0=">TaB|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,10,1)=="D"){$0=">TaD|"substr($1,2,length($1))}} 1' > Triticum_aestivum.IWGSC.pep.allname2.fa
###CS /data1/home/xuebo/wheat/Ref/pep
seqkit sort -n Triticum_aestivum.IWGSC.pep.all.fa.gz > Triticum_aestivum.IWGSC.pep.allsort.fa
grep -n "TraesCSU" Triticum_aestivum.IWGSC.pep.allsort.fa | head -n 1  1161346
sed -n '1,1161345p' Triticum_aestivum.IWGSC.pep.allsort.fa > Triticum_aestivum.IWGSC.pep.allsort2.fa
cat Triticum_aestivum.IWGSC.pep.allsort2.fa | awk '/>/ {$0 =substr($1,1,19)} 1' > Triticum_aestivum.IWGSC.pep.allname.fa
cat Triticum_aestivum.IWGSC.pep.allname.fa | awk '/>/ {if(substr($1,10,1)=="A"){$0=">TaA|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,10,1)=="B"){$0=">TaB|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,10,1)=="D"){$0=">TaD|"substr($1,2,length($1))}} 1' > CS_ABD.pep.fa
###urartu /data1/home/xuebo/wheat/othergenome/urartu/pep
seqkit sort -n Triticum_urartu.ASM34745v1.pep.all.fa.gz > Triticum_urartu.ASM34745v1.pep.allsort.fa
cat Triticum_urartu.ASM34745v1.pep.allsort.fa | awk '/>/ {$0 =substr($1,1,length($1))} 1' > Triticum_urartu.ASM34745v1.pep.allname.fa
cat Triticum_urartu.ASM34745v1.pep.allname.fa | awk '/>/ {$0=">Tu|"substr($1,2,length($1))} 1'  > Tu.pep.fa
###barley /data1/home/xuebo/wheat/othergenome/barley/pep
seqkit sort -n Hordeum_vulgare.IBSC_v2.pep.all.fa.gz > Hordeum_vulgare.IBSC_v2.pep.allsort.fa
cat Hordeum_vulgare.IBSC_v2.pep.allsort.fa | awk '/>/ {$0 =substr($1,1,length($1))} 1' > Hordeum_vulgare.IBSC_v2.pep.allname.fa
cat Hordeum_vulgare.IBSC_v2.pep.allname.fa | awk '/>/ {$0=">Hv|"substr($1,2,length($1))} 1'  > Hv.pep.fa
###tauschii /data1/home/xuebo/wheat/othergenome/Aetgenome/pep
seqkit sort -n Aegilops_tauschii.Aet_v4.0.pep.all.fa.gz > Aegilops_tauschii.Aet_v4.0.pep.allsort.fa
cat Aegilops_tauschii.Aet_v4.0.pep.allsort.fa | awk '/>/ {$0 =substr($1,1,length($1))} 1' > Aegilops_tauschii.Aet_v4.0.pep.allname.fa
cat Aegilops_tauschii.Aet_v4.0.pep.allname.fa | awk '/>/ {$0=">Ae|"substr($1,2,length($1))} 1'  > Ae.pep.fa
###wildemmer  /data1/home/xuebo/wheat/othergenome/dicoccoides/pep
seqkit sort -n Triticum_dicoccoides.WEWSeq_v.1.0.pep.all.fa.gz > Triticum_dicoccoides.WEWSeq_v.1.0.pep.allsort.fa
cat Triticum_dicoccoides.WEWSeq_v.1.0.pep.allsort.fa | awk '/>/ {$0=substr($1,1,length($1))} 1' > Triticum_dicoccoides.WEWSeq_v.1.0.pep.allname.fa 
cat Triticum_dicoccoides.WEWSeq_v.1.0.pep.allname.fa | awk '/>/ {if(substr($1,8,1)=="A"){$0=">TdA|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,8,1)=="B"){$0=">TdB|"substr($1,2,length($1))}} 1'  > wildemmer_AB.pep.fa
###Durum /data1/home/xuebo/wheat/othergenome/Durum/pep
seqkit sort -n Triticum_turgidum.Svevo.v1.pep.all.fa.gz > Triticum_turgidum.Svevo.v1.pep.allsort.fa
grep -n "TRITD1A"  Triticum_turgidum.Svevo.v1.pep.allsort.fa | head -n 1
sed -n '47954,1797930p' Triticum_turgidum.Svevo.v1.pep.allsort.fa > Triticum_turgidum.Svevo.v1.pep.allsort2.fa
cat Triticum_turgidum.Svevo.v1.pep.allsort2.fa | awk '/>/ {$0=substr($1,1,length($1))} 1' > Triticum_turgidum.Svevo.v1.pep.allname.fa
cat Triticum_turgidum.Svevo.v1.pep.allname.fa | awk '/>/ {if(substr($1,8,1)=="A"){$0=">TtA|"substr($1,2,length($1))}} 1' | awk '/>/ {if(substr($1,8,1)=="B"){$0=">TtB|"substr($1,2,length($1))}} 1'  > Durum_AB.pep.fa
##合起来 /data1/home/xuebo/Projects/Speciation/species_split_time/SwiftOrtho
cat /data1/home/xuebo/wheat/Ref/pep/CS_ABD.pep.fa  /data1/home/xuebo/wheat/othergenome/urartu/pep/Tu.pep.fa /data1/home/xuebo/wheat/othergenome/barley/pep/Hv.pep.fa /data1/home/xuebo/wheat/othergenome/Aetgenome/pep/Ae.pep.fa  /data1/home/xuebo/wheat/othergenome/dicoccoides/pep/wildemmer_AB.pep.fa  /data1/home/xuebo/wheat/othergenome/Durum/pep/Durum_AB.pep.fa > species_all5.pep.fa
nohup python /data1/home/xuebo/software/SwiftOrtho/bin/find_hit.py -p blastp -i species_all5.pep.fa -d species_all5.pep.fa -o species_all5.pep.fa.sc -e 1e-5 -s 111111 -a 120 > nohupout 2>& 1 &
####在204上面用不同的软件做
diamond makedb --in species_all5.pep.fa -d nr
nohup diamond blastp -d nr -q species_all5.pep.fa -o matches.m8 -p 120 > nohupout 2>& 1 &
nohup python /data1/home/xuebo/software/SwiftOrtho/bin/find_orth.py -i matches.m8 -c 0.5 -y 0 -a 120 > input.fsa.sc.orth > nohupout2 2>& 1 &
nohup python /data1/home/xuebo/software/SwiftOrtho/bin/find_cluster.py -i input.fsa.sc.orth -a mcl -I 1.5 > input.fsa.sc.orth.mcl  | grep -n "Hv"> nohupout3 2>& 1 &
nohup python /data1/home/xuebo/software/SwiftOrtho/bin/find_orth.py -i matches.m8 -c 0.5 -y 50 -a 120 > input.fsa.sc.orth2 > nohupout2 2>& 1 &
nohup python /data1/home/xuebo/software/SwiftOrtho/bin/find_cluster.py -i input.fsa.sc.orth -a mcl -I 1.5 > input.fsa.sc.orth.mcl2 > nohupout3 2>& 1 &
#####
cat input.fsa.sc.orth.mcl  | grep -n "Hv" && grep -n "Ae" && grep -n "Tu" && grep -n "TdA" && grep -n "TdB" && grep -n "TtA" && grep -n "TtB" && grep -n "TaA" && grep -n "TaB" && grep -n "TaD"
cat input.fsa.sc.orth.mcl  | grep -n "Hv" | grep -n "Ae" | grep -n "Tu" | grep -n "TdA" | grep -n "TdB" | grep -n "TtA" | grep -n "TtB" | grep -n "TaA" | grep -n "TaB" | grep -n "TaD"  > singlecopygene.txt
cat input.fsa.sc.orth.mcl  | grep "Hv" | grep "Ae" | grep "Tu" | grep "TdA" | grep "TdB" | grep "TtA" | grep "TtB" | grep "TaA" | grep "TaB" | grep "TaD"  > singlecopygene.txt
###现在是挑选个体，每个亚种挑选5个 /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta
cp a_iwgscV1_forscan.fa.gz a_iwgscV1_forGATK.fa.gz
samtools faidx a_iwgscV1_forGATK.fa.gz     				bgzip才可以
/data1/home/xuebo/software/gatk-4.1.4.0/gatk CreateSequenceDictionary -R a_iwgscV1_forGATK.fa.gz 
cp d_iwgscV1_forscan.fa.gz d_iwgscV1_forGATK.fa.gz   
samtools faidx d_iwgscV1_forGATK.fa.gz     				
/data1/home/xuebo/software/gatk-4.1.4.0/gatk CreateSequenceDictionary -R d_iwgscV1_forGATK.fa.gz 
先把barley的bam文件和一下 /data2/xuebo/Projects/Speciation/tree/barleybam
samtools merge -h chr1.sorted.bam -c -p -@ 10 barley_Alineage.bam  chr1.sorted.bam chr2.sorted.bam chr7.sorted.bam chr8.sorted.bam chr13.sorted.bam chr14.sorted.bam chr19.sorted.bam chr20.sorted.bam chr25.sorted.bam chr26.sorted.bam chr31.sorted.bam chr32.sorted.bam chr37.sorted.bam chr38.sorted.bam &
samtools merge -h chr3.sorted.bam -c -p -@ 10 barley_Blineage.bam  chr3.sorted.bam chr4.sorted.bam chr9.sorted.bam chr10.sorted.bam chr15.sorted.bam chr16.sorted.bam chr21.sorted.bam chr22.sorted.bam chr27.sorted.bam chr28.sorted.bam chr33.sorted.bam chr34.sorted.bam chr39.sorted.bam chr40.sorted.bam &
samtools merge -h chr5.sorted.bam -c -p -@ 10 barley_Dlineage.bam  chr5.sorted.bam chr6.sorted.bam chr11.sorted.bam chr12.sorted.bam chr17.sorted.bam chr18.sorted.bam chr23.sorted.bam chr24.sorted.bam chr29.sorted.bam chr30.sorted.bam chr35.sorted.bam chr36.sorted.bam chr41.sorted.bam chr42.sorted.bam &
***barleyA  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/barleyA
#!/bin/bash 
for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data2/xuebo/Projects/Speciation/tree/barleybam/barley_Alineage.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O barley_A_${i}.g.vcf 
	  cat barley_A_${i}.g.vcf  | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > barley_A_${i}.fa 
done
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***wild_einkorn  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/wild_einkorn     
#!/bin/bash 
for taxa in {"A_0001","A_0003","A_0005","A_0012","A_0028"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/A/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_A_${i}.g.vcf 
	  cat ${taxa}_A_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_A_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
#!/bin/bash
for taxa in {"A_0001","A_0003","A_0005","A_0012","A_0028"}
do
	cd ${taxa}
	for ii in *.fa; 
	do
			ii=${ii%.fa};
			sed -i -e '$a\' ${ii}.fa
	done
	cd ..
done
***'urartu  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/urartu     
#!/bin/bash 
for taxa in {"A_0066","A_0070","A_0071","A_0076","A_0085"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/A/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_A_${i}.g.vcf 
	  cat ${taxa}_A_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_A_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***wild_emmerA  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/wild_emmerA     
#!/bin/bash 
for taxa in {"AB_0023","AB_0024","AB_0034","AB_0037","AB_0045"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/AB/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_A_${i}.g.vcf 
	  cat ${taxa}_A_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_A_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***durumA  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/durumA     
#!/bin/bash 
for taxa in {"AB_0114","AB_0115","AB_0116","AB_0119","AB_0122"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/AB/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_A_${i}.g.vcf 
	  cat ${taxa}_A_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_A_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***landraceA  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/landraceA     
#!/bin/bash 
for taxa in {"ABD_0066","ABD_0069","ABD_0077","ABD_0113","ABD_0165"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/a_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/ABD/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100A/singlecopygene_A_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_A_${i}.g.vcf 
	  cat ${taxa}_A_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_A_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***barleyB  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/barleyB
#!/bin/bash 
for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz \
	  -I /data2/xuebo/Projects/Speciation/tree/barleybam/barley_Blineage.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100B/singlecopygene_B_${i}.bed -ERC BP_RESOLUTION -O barley_B_${i}.g.vcf 
	  cat barley_B_${i}.g.vcf  | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > barley_B_${i}.fa 
done
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***speltoides  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/speltoides     
#!/bin/bash 
for taxa in {"S002","S003","S006","S007","S009"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz \
	  -I /data2/xuebo/Projects/Speciation/bamS10/ByChrMergeSforHaplotypeCaller/${taxa}.rmdup.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100B/singlecopygene_B_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_B_${i}.g.vcf 
	  cat ${taxa}_B_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_B_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***wild_emmerB  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/wild_emmerB     
#!/bin/bash 
for taxa in {"AB_0023","AB_0024","AB_0034","AB_0037","AB_0045"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/AB/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100B/singlecopygene_B_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_B_${i}.g.vcf 
	  cat ${taxa}_B_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_B_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***durumB  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/durumB     
#!/bin/bash 
for taxa in {"AB_0114","AB_0115","AB_0116","AB_0116","AB_0122"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/AB/${taxa}.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100B/singlecopygene_B_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_B_${i}.g.vcf 
	  cat ${taxa}_B_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_B_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***landraceB  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/landraceB
#!/bin/bash 
for taxa in {"ABD_0066","ABD_0069","ABD_0077","ABD_0113","ABD_0165"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/ABD/${taxa}.bam  --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100B/singlecopygene_B_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_B_${i}.g.vcf 
	  cat ${taxa}_B_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_B_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***barleyD  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/barleyD
#!/bin/bash 
for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/d_iwgscV1_forGATK.fa.gz \
	  -I /data2/xuebo/Projects/Speciation/tree/barleybam/barley_Dlineage.bam --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100D/singlecopygene_D_${i}.bed -ERC BP_RESOLUTION -O barley_D_${i}.g.vcf 
	  cat barley_D_${i}.g.vcf  | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > barley_D_${i}.fa 
done
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***tauschii  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/tauschii
#!/bin/bash 
for taxa in {"D_0003","D_0004","D_0008","D_0021","D_0025"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/d_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/D/${taxa}.bam  --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100D/singlecopygene_D_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_D_${i}.g.vcf 
	  cat ${taxa}_D_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_D_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***strangulata  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/strangulata
#!/bin/bash 
for taxa in {"D_0001","D_0010","D_0014","D_0015","D_0028"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/d_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/D/${taxa}.bam  --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100D/singlecopygene_D_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_D_${i}.g.vcf 
	  cat ${taxa}_D_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_D_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
***landraceD  getgatkcaller.sh    /data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/landraceD
#!/bin/bash 
for taxa in {"ABD_0066","ABD_0069","ABD_0077","ABD_0113","ABD_0165"}
do
	echo "$taxa is running..." 
	mkdir ${taxa}
	cd ${taxa}
	for i in {1..100}
	do
	  /data1/home/xuebo/software/gatk-4.1.4.0/gatk --java-options "-Xmx35g" HaplotypeCaller  -R /data1/home/xuebo/Project/reference/d_iwgscV1_forGATK.fa.gz \
	  -I /data3/wgs/bam/ABD/${taxa}.bam  --intervals /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1.gene100D/singlecopygene_D_${i}.bed -ERC BP_RESOLUTION -O ${taxa}_D_${i}.g.vcf 
	  cat ${taxa}_D_${i}.g.vcf | perl /data1/home/xuebo/software/haploblocks/tools/gvcf2fasta_nogaps.pl > ${taxa}_D_${i}.fa 
  done
  cd ..
done  
nohup sh getgatkcaller.sh > nohupgatkcaller  2>& 1 & 
######################################################现在想要做A的情况
##首先做的是文件合并 /data2/xuebo/Projects/Speciation/species_split_time/alin_fasta/A/fasta_file  cat.sh
#!/bin/bash
infile="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta"
inwild_einkorn="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/wild_einkorn"
inurartu="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/urartu"
inwild_emmerA="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/wild_emmerA"
indurumA="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/durumA"
inlandraceA="/data2/xuebo/Projects/Speciation/species_split_time/gvcf2fasta/landraceA"
for i in {1..100}
do
	cat ${infile}/barleyA/barley_A_${i}.fa   ${inwild_einkorn}/A_0001/A_0001_A_${i}.fa  ${inwild_einkorn}/A_0003/A_0003_A_${i}.fa  ${inwild_einkorn}/A_0005/A_0005_A_${i}.fa ${inwild_einkorn}/A_0012/A_0012_A_${i}.fa  ${inwild_einkorn}/A_0028/A_0028_A_${i}.fa \
	${inurartu}/A_0066/A_0066_A_${i}.fa  ${inurartu}/A_0070/A_0070_A_${i}.fa  ${inurartu}/A_0071/A_0071_A_${i}.fa ${inurartu}/A_0076/A_0076_A_${i}.fa  ${inurartu}/A_0085/A_0085_A_${i}.fa \
	${inwild_emmerA}/AB_0023/AB_0023_A_${i}.fa  ${inwild_emmerA}/AB_0024/AB_0024_A_${i}.fa  ${inwild_emmerA}/AB_0034/AB_0034_A_${i}.fa ${inwild_emmerA}/AB_0037/AB_0037_A_${i}.fa  ${inwild_emmerA}/AB_0045/AB_0045_A_${i}.fa \
	${indurumA}/AB_0114/AB_0114_A_${i}.fa  ${indurumA}/AB_0115/AB_0115_A_${i}.fa  ${indurumA}/AB_0116/AB_0116_A_${i}.fa ${indurumA}/AB_0119/AB_0119_A_${i}.fa  ${indurumA}/AB_0122/AB_0122_A_${i}.fa \
	${inlandraceA}/ABD_0066/ABD_0066_A_${i}.fa  ${inlandraceA}/ABD_0069/ABD_0069_A_${i}.fa  ${inlandraceA}/ABD_0077/ABD_0077_A_${i}.fa ${inlandraceA}/ABD_0113/ABD_0113_A_${i}.fa  ${inlandraceA}/ABD_0165/ABD_0165_A_${i}.fa > A_${i}.alin.fa
done 
###########beast
barley
A001
A003   
A005
A012
A028
A066
A070
A071
A076
A085
B023
B024
B034
B037
B045
B114
B115
B116
B119
B122
FIN-L1
GEO-L1
IRN-L3
UZB-L1
HCM
***产生.xml文件
#!/bin/bash
for i in {1..100}
do
  java -jar /data2/xuebo/Projects/Speciation/javaCode/C39_get_xmlfile_forbeast.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/alin_fasta/A/fasta_file/A_${i}.alin.fa \
  --file2 /data2/xuebo/Projects/Speciation/species_split_time/alin_fasta/A/fasta_file/A_2.xml --out A_${i}.xml 
done
#!/bin/bash
for i in {1..100}
do
	/data1/home/xuebo/software/beast/bin/beast -threads 120 /data2/xuebo/Projects/Speciation/species_split_time/alin_fasta/A/xml_file/A_${i}.xml &
	monitor java 24 10s
done
#!/bin/bash
for i in {1..100}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 50 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/alin_fasta/A/runbeast/A_${i}.trees A_${i}.mcc.tre &
	monitor java 24 2s
done
#####改名字
#!/bin/bash
for i in {1..100}
do
    echo "TREE"${i};
		sed -n '63,63p' A_${i}.mcc.tre  |  sed 's/TREE1/'"TREE${i}"'/' >> Alineage_1_100.mcc.tre
done

for i in {1..100};do echo "TREE"${i};sed -n '63,63p' A_${i}.mcc.tre  |  sed 's/TREE1/TREE${i}/' >> Alineage_1_100.mcc.tre;done

######################################################以上流程跑通啦，现在要做的是对assembly的基因进行分析，算时间
#######准备文件
****现在是有一个1：1而且>1000bp的基因的list,找到assembly文件里面ABD都有的基因，之后复制到指定路径
#!/bin/bash
java -jar /data2/xuebo/Projects/Speciation/javaCode/C42_move1_1genefile.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1_more1000.txt --file2 /data1/home/xuebo/SRAssembler/assembly --out /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/singlecopygene1_1_more1000_assembly.txt &
monitor java 22 2s
*****现在得到了singlecopygene1_1_more1000_assembly.txt，把barley的基因找出来，得到bed文件，之后
*https://bioweb.pasteur.fr/docs/modules/bedtools/2.17.0/content/fastafromBed.html
* fastaFromBed -fi Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa -bed test.bed -fo test.fa.out -name
* chr1H	272379  275627	HORVU1Hr1G000110
****Hvbed.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$5}' > nameA.txt
    cat test111.txt    | awk '{print$6}' > nameB.txt
    cat test111.txt    | awk '{print$7}' > nameD.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
    nameB=$(sed -n 1p nameB.txt)
    echo $nameB
    nameD=$(sed -n 1p nameD.txt)
    echo $nameD
    cat test111.txt | awk '{printf $1 "\t" $2 "\t" $3 "\tbarley\n"}'  >  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A/$nameA/temp_barley.bed
    cat test111.txt | awk '{printf $1 "\t" $2 "\t" $3 "\tbarley\n"}'  >  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B/$nameB/temp_barley.bed
    cat test111.txt | awk '{printf $1 "\t" $2 "\t" $3 "\tbarley\n"}'  >  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D/$nameD/temp_barley.bed
done
****getHvfaA.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$5}' > nameA.txt
    #cat test111.txt    | awk '{print$6}' > nameB.txt
    #cat test111.txt    | awk '{print$7}' > nameD.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A/$nameA
    		rm barley_${nameA}.fa
				fastaFromBed -fi /data1/home/xuebo/Project/reference/orthergenome/barley/Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa -bed temp_barley.bed -fo temp_barley.fa -name
    		awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' temp_barley.fa > barley_${nameA}.fasta
    cd ../../
done
****getHvfaB.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
    cat test111.txt    | awk '{print$6}' > nameB.txt
		nameB=$(sed -n 1p nameB.txt)
    echo $nameB
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B/$nameB
				fastaFromBed -fi /data1/home/xuebo/Project/reference/orthergenome/barley/Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa -bed temp_barley.bed -fo temp_barley.fa -name
    		awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' temp_barley.fa > barley_${nameB}.fasta
    cd ../../
done
****getHvfaD.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
    cat test111.txt    | awk '{print$7}' > nameD.txt
		nameD=$(sed -n 1p nameD.txt)
    echo $nameD
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D/$nameD
				fastaFromBed -fi /data1/home/xuebo/Project/reference/orthergenome/barley/Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa -bed temp_barley.bed -fo temp_barley.fa -name
    		awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' temp_barley.fa > barley_${nameD}.fasta
    cd ../../
done
#######多序列比对，使用muscle是，生成多序列比对文件
**首先要做的是改名字
****renameA.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$5}' > nameA.txt
    #cat test111.txt    | awk '{print$6}' > nameB.txt
    #cat test111.txt    | awk '{print$7}' > nameD.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A/$nameA
    		sed -e 's/Landrace/LandraceA/' Landrace_${nameA}.fasta > A_Landrace_${nameA}.fasta
    		sed -e 's/Club_wheat/ClubA/' Club_wheat_${nameA}.fasta > A_Club_wheat_${nameA}.fasta
    		sed -e 's/Cultivar/CultivarA/' Cultivar_${nameA}.fasta > A_Cultivar_${nameA}.fasta
    		sed -e 's/Domesticated_emmer/DomemmerA/' Domesticated_emmer_${nameA}.fasta > A_Domesticated_emmer_${nameA}.fasta
    		sed -e 's/Durum/DurumA/' Durum_${nameA}.fasta > A_Durum_${nameA}.fasta
    		sed -e 's/Georgian_wheat/GeorgianA/' Georgian_wheat_${nameA}.fasta > A_Georgian_wheat_${nameA}.fasta
    		sed -e 's/Indian_dwarf_wheat/IndianA/' Indian_dwarf_wheat_${nameA}.fasta > A_Indian_dwarf_wheat_${nameA}.fasta
    		sed -e 's/Ispahanicum/IspahanA/' Ispahanicum_${nameA}.fasta > A_Ispahanicum_${nameA}.fasta
    		sed -e 's/Khorasan_wheat/KhorasanA/' Khorasan_wheat_${nameA}.fasta > A_Khorasan_wheat_${nameA}.fasta
    		sed -e 's/Macha/MachaA/' Macha_${nameA}.fasta > A_Macha_${nameA}.fasta
    		sed -e 's/Persian_wheat/PersianA/' Persian_wheat_${nameA}.fasta > A_Persian_wheat_${nameA}.fasta
    		sed -e 's/Polish_wheat/PolishA/' Polish_wheat_${nameA}.fasta > A_Polish_wheat_${nameA}.fasta
    		sed -e 's/Rivet_wheat/RivetA/' Rivet_wheat_${nameA}.fasta > A_Rivet_wheat_${nameA}.fasta
    		sed -e 's/Spelt/SpeltA/' Spelt_${nameA}.fasta > A_Spelt_${nameA}.fasta
    		#sed -e 's/Synthetic/SyntheticA/' Synthetic_${nameA}.fasta > A_Synthetic_${nameA}.fasta
    		sed -e 's/Tibetan_semi_wild/TibetanA/' Tibetan_semi_wild_${nameA}.fasta > A_Tibetan_semi_wild_${nameA}.fasta
    		sed -e 's/Vavilovii/VaviloviiA/' Vavilovii_${nameA}.fasta > A_Vavilovii_${nameA}.fasta
    		sed -e 's/Wild_emmer/WildemmerA/' Wild_emmer_${nameA}.fasta > A_Wild_emmer_${nameA}.fasta
    		sed -e 's/Xinjiang_wheat/XinjiangA/' Xinjiang_wheat_${nameA}.fasta > A_Xinjiang_wheat_${nameA}.fasta
    		sed -e 's/Yunan_wheat/YunanA/' Yunan_wheat_${nameA}.fasta > A_Yunan_wheat_${nameA}.fasta
    		sed -e 's/Domesticated_einkorn/Domeink/' Domesticated_einkorn_${nameA}.fasta > A_Domesticated_einkorn_${nameA}.fasta
    		sed -e 's/Urartu/Urartu/' Urartu_${nameA}.fasta > A_Urartu_${nameA}.fasta
    		sed -e 's/Wild_einkorn/Wildeink/' Wild_einkorn_${nameA}.fasta > A_Wild_einkorn_${nameA}.fasta
    		cat A_*fasta > subA_${nameA}.fasta
    	cd ../../
done
****renameB.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
    cat test111.txt    | awk '{print$6}' > nameB.txt
		nameB=$(sed -n 1p nameB.txt)
    echo $nameB
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B/$nameB
    		sed -e 's/Landrace/LandraceB/' Landrace_${nameB}.fasta > B_Landrace_${nameB}.fasta
    		sed -e 's/Club_wheat/ClubB/' Club_wheat_${nameB}.fasta > B_Club_wheat_${nameB}.fasta
    		sed -e 's/Cultivar/CultivarB/' Cultivar_${nameB}.fasta > B_Cultivar_${nameB}.fasta
    		sed -e 's/Domesticated_emmer/DomemmerB/' Domesticated_emmer_${nameB}.fasta > B_Domesticated_emmer_${nameB}.fasta
    		sed -e 's/Durum/DurumB/' Durum_${nameB}.fasta > B_Durum_${nameB}.fasta
    		sed -e 's/Georgian_wheat/GeorgianB/' Georgian_wheat_${nameB}.fasta > B_Georgian_wheat_${nameB}.fasta
    		sed -e 's/Indian_dwarf_wheat/IndianB/' Indian_dwarf_wheat_${nameB}.fasta > B_Indian_dwarf_wheat_${nameB}.fasta
    		sed -e 's/Ispahanicum/IspahanB/' Ispahanicum_${nameB}.fasta > B_Ispahanicum_${nameB}.fasta
    		sed -e 's/Khorasan_wheat/KhorasanB/' Khorasan_wheat_${nameB}.fasta > B_Khorasan_wheat_${nameB}.fasta
    		sed -e 's/Macha/MachaB/' Macha_${nameB}.fasta > B_Macha_${nameB}.fasta
    		sed -e 's/Persian_wheat/PersianB/' Persian_wheat_${nameB}.fasta > B_Persian_wheat_${nameB}.fasta
    		sed -e 's/Polish_wheat/PolishB/' Polish_wheat_${nameB}.fasta > B_Polish_wheat_${nameB}.fasta
    		sed -e 's/Rivet_wheat/RivetB/' Rivet_wheat_${nameB}.fasta > B_Rivet_wheat_${nameB}.fasta
    		sed -e 's/Spelt/SpeltB/' Spelt_${nameB}.fasta > B_Spelt_${nameB}.fasta
    		#sed -e 's/Synthetic/SyntheticB/' Synthetic_${nameB}.fasta > B_Synthetic_${nameB}.fasta
    		sed -e 's/Tibetan_semi_wild/TibetanB/' Tibetan_semi_wild_${nameB}.fasta > B_Tibetan_semi_wild_${nameB}.fasta
    		sed -e 's/Vavilovii/VaviloviiB/' Vavilovii_${nameB}.fasta > B_Vavilovii_${nameB}.fasta
    		sed -e 's/Wild_emmer/WildemmerB/' Wild_emmer_${nameB}.fasta > B_Wild_emmer_${nameB}.fasta
    		sed -e 's/Xinjiang_wheat/XinjiangB/' Xinjiang_wheat_${nameB}.fasta > B_Xinjiang_wheat_${nameB}.fasta
    		sed -e 's/Yunan_wheat/YunanB/' Yunan_wheat_${nameB}.fasta > B_Yunan_wheat_${nameB}.fasta
    		sed -e 's/Speltoides/SpeltoidB/' Speltoides_${nameB}.fasta > B_Speltoides_${nameB}.fasta
    		cat B_*fasta > subB_${nameB}.fasta
    	cd ../../
done
****renameD.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > test111.txt     
    cat test111.txt    | awk '{print$7}' > nameD.txt
		nameD=$(sed -n 1p nameD.txt)
    echo $nameD
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D/$nameD
    		sed -e 's/Landrace/LandraceD/' Landrace_${nameD}.fasta > D_Landrace_${nameD}.fasta
    		sed -e 's/Club_wheat/ClubD/' Club_wheat_${nameD}.fasta > D_Club_wheat_${nameD}.fasta
    		sed -e 's/Cultivar/CultivarD/' Cultivar_${nameD}.fasta > D_Cultivar_${nameD}.fasta
    		sed -e 's/Indian_dwarf_wheat/IndianD/' Indian_dwarf_wheat_${nameD}.fasta > D_Indian_dwarf_wheat_${nameD}.fasta
    		sed -e 's/Macha/MachaD/' Macha_${nameD}.fasta > D_Macha_${nameD}.fasta
    		sed -e 's/Spelt/SpeltD/' Spelt_${nameD}.fasta > D_Spelt_${nameD}.fasta
    		#sed -e 's/Synthetic/SyntheticD/' Synthetic_${nameD}.fasta > D_Synthetic_${nameD}.fasta
    		sed -e 's/Tibetan_semi_wild/TibetanD/' Tibetan_semi_wild_${nameD}.fasta > D_Tibetan_semi_wild_${nameD}.fasta
    		sed -e 's/Vavilovii/VaviloviiD/' Vavilovii_${nameD}.fasta > D_Vavilovii_${nameD}.fasta
    		sed -e 's/Xinjiang_wheat/XinjiangD/' Xinjiang_wheat_${nameD}.fasta > D_Xinjiang_wheat_${nameD}.fasta
    		sed -e 's/Yunan_wheat/YunanD/' Yunan_wheat_${nameD}.fasta > D_Yunan_wheat_${nameD}.fasta
    		sed -e 's/Anathera/Anathera/' Anathera_${nameD}.fasta > D_Anathera_${nameD}.fasta
    		sed -e 's/Meyeri/Meyeri/' Meyeri_${nameD}.fasta > D_Meyeri_${nameD}.fasta
    		sed -e 's/Strangulata/Strang/' Strangulata_${nameD}.fasta > D_Strangulata_${nameD}.fasta
    		cat D_*fasta > subD_${nameD}.fasta
    	cd ../../
done
###比对
****muscleA.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > testmuscle.txt     
		cat testmuscle.txt   | awk '{print$5}' > nameA.txt
    #cat test111.txt    | awk '{print$6}' > nameB.txt
    #cat test111.txt    | awk '{print$7}' > nameD.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A/$nameA
				muscle -in subA_${nameA}.fasta -out /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A_alin/subA_${nameA}.afa -maxiters 2 &
				monitor muscle 100 5s
    cd ../../
done
****muscleB.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > testmuscleB.txt     
		cat testmuscleB.txt   | awk '{print$6}' > nameB.txt
		nameB=$(sed -n 1p nameB.txt)
    echo $nameB
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B/$nameB
				muscle -in subB_${nameB}.fasta -out /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B_alin/subB_${nameB}.afa -maxiters 2 &
				monitor muscle 60 5s
    cd ../../
done
****muscleD.sh
#!/bin/bash
cat singlecopygene1_1_more1000_assembly_addHv.txt | while read line
do
    echo $line > testmuscleD.txt     
		cat testmuscleD.txt   | awk '{print$7}' > nameD.txt
		nameD=$(sed -n 1p nameD.txt)
    echo $nameD
    cd /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D/$nameD
				muscle -in subD_${nameD}.fasta -out /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D_alin/subD_${nameD}.afa -maxiters 2 &
				monitor muscle 160 5s
    cd ../../
done
####B_alin里面的文件个数是4975，D_alin里面的文件是4976，所以要检查一下哪个文件没有
find /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D_alin -name *.afa > /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/listB.txt
TraesCS1D02G221000 这个里面没有组装的结果
chr1H	410438725	410442190	HORVU1Hr1G056020	TraesCS1A02G219400	TraesCS1B02G232800	TraesCS1D02G221000 这一行是不能用的
TraesCS3B02G245300
TraesCS3B02G315800
chr3H	381156059	381160615	HORVU3Hr1G052470	TraesCS3A02G215300	TraesCS3B02G245300	TraesCS3D02G217300
chr3H	525392337	525395803	HORVU3Hr1G069290	TraesCS3A02G282100	TraesCS3B02G315800	TraesCS3D02G282100
所以A的文件夹里面要删除
TraesCS1A02G219400
TraesCS3A02G215300
TraesCS3A02G282100
所以B的文件夹里面要删除
TraesCS1B02G232800
TraesCS3B02G245300
TraesCS3B02G315800
所以D的文件夹里面要删除
TraesCS1D02G221000
TraesCS3D02G217300
TraesCS3D02G282100
改好的文件夹在这里/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/singlecopygene1_1_more1000_assembly_addHv2_4974.txt
####现在做的是物种的比对
****catspeciesfa.sh
#!/bin/bash
i=0
cat singlecopygene1_1_more1000_assembly_addHv2_4974.txt | while read line
do
    i=$((i+1))
    echo $i
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$5}' > nameA.txt
    cat test111.txt    | awk '{print$6}' > nameB.txt
    cat test111.txt    | awk '{print$7}' > nameD.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
		nameB=$(sed -n 1p nameB.txt)
    echo $nameB
		nameD=$(sed -n 1p nameD.txt)
    echo $nameD 
		cat ./B/${nameB}/barley_${nameB}.fasta  ./A/${nameA}/A_Wild_einkorn_${nameA}.fasta  ./A/${nameA}/A_Urartu_${nameA}.fasta  ./A/${nameA}/A_Wild_emmer_${nameA}.fasta  ./A/${nameA}/A_Rivet_wheat_${nameA}.fasta  ./A/${nameA}/A_Landrace_${nameA}.fasta  \
		./B/${nameB}/B_Speltoides_${nameB}.fasta  ./B/${nameB}/B_Wild_emmer_${nameB}.fasta  ./B/${nameB}/B_Rivet_wheat_${nameB}.fasta  ./B/${nameB}/B_Landrace_${nameB}.fasta  ./D/${nameD}/D_Meyeri_${nameD}.fasta  ./D/${nameD}/D_Strangulata_${nameD}.fasta  ./D/${nameD}/D_Landrace_${nameD}.fasta > ./all_alin/species_$i.fasta
done
****musclespecies.sh
#!/bin/bash
for i in {1..4974}
do 
	muscle -in ./all_alin/species_$i.fasta -out ./all_alin/species_$i.afa -maxiters 2 &
  monitor muscle 120 5s
done
find /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/all_alin/alin -name *.afa > /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/listall.txt
nohup muscle -in species_67.fasta -out ./alin/species_67.afa -maxiters 2  > nohup67 2>& 1 &
nohup muscle -in species_391.fasta -out ./alin/species_391.afa -maxiters 2  > nohup391 2>& 1 &
nohup muscle -in species_3862.fasta -out ./alin/species_3862.afa -maxiters 2  > nohup3862 2>& 1 &
67 chr1H	332740736	332751652	HORVU1Hr1G045720	TraesCS1A02G174000	TraesCS1B02G190400	TraesCS1D02G164800
391 chr2H	716583162	716885241	HORVU2Hr1G109430	TraesCS2A02G490000	TraesCS2B02G518300	TraesCS2D02G490400
3862 chr3H	31954734	31955752	HORVU3Hr1G014090	TraesCS3A02G077900	TraesCS3B02G092800	TraesCS3D02G078500
#!/bin/bash
for i in {68..390}
do 
	j=$((i-1))
	cp ./alin/species_$i.afa  ./alin2/species_$j.afa 
done
#!/bin/bash
for i in {392..3861}
do 
	j=$((i-2))
	cp ./alin/species_$i.afa  ./alin2/species_$j.afa 
done
#!/bin/bash
for i in {3863..4974}
do 
	j=$((i-3))
	cp ./alin/species_$i.afa  ./alin2/species_$j.afa 
done
**改成一行的  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < species_9.afa
cat species_9.afa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > species1line_9.alin.fa
sed -i '1d' species1line_9.alin.fa
#!/bin/bash
for i in {1..4971}
do 
	cat species_$i.afa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > species1line_$i.alin.fa
	sed -i '1d' species1line_$i.alin.fa
done
***产生.xml文件
#!/bin/bash
for i in {1..4971}
do
  java -jar /data2/xuebo/Projects/Speciation/javaCode/C39_get_xmlfile_forbeast2.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/all_alin/alin2/species1line_$i.alin.fa \
  --file2 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/species1line_6.xml --out ./xml_species/species1line_${i}.xml
done
********runbeast.sh
#!/bin/bash
for i in {1..4971}
do
	/data1/home/xuebo/software/beast/bin/beast -threads 10 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/xml_species/species1line_${i}.xml &
	monitor java 22 10s
done
********runtreeannotator.sh
/data1/home/xuebo/software/beast/bin/treeannotator -burnin 50 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_species/species1line_6.trees species1line_6.mcc.tre 
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 20 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_species/species1line_${i}.trees species1line_${i}.mcc.tre &
	monitor java 24 2s
done
*****改名字
#!/bin/bash
for i in {1..4971}
do
    echo "TREE"${i};
		sed -n '37,37p' species1line_${i}.mcc.tre |  sed 's/TREE1/'"TREE${i}"'/' >>  species1line1_4971.mcc.tre
done
*******runAPE.r
#!/usr/bin/Rscript
library(ape)
setwd("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_species")
DD = read.table("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/treeannotator_species/species1line1_4971.mcc.tre",header = F, sep = "\t",quote = "\"")
allapetree <- matrix("NA", nrow = 4971, ncol = 1)
for(i in 1:4971){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol =  6.9e-3,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  write.tree(res,paste("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_species/apeTree/",i,".tree",sep=""))
}
##最后测试了6次，tol =  1e-3,这个结果最为合理可靠，最后使用这个画图整理
cat *tree > ../ape6_species1line1_4971.mcc.tre
shuf -n 1500 ape6_species1line1_4971.mcc.tre > ape6_shuf1500_species1line1_4971.mcc.tre
**算时间
/data1/home/xuebo/software/beast/bin/treeannotator -burnin 10 -heights median  species1line1_4971.mcc.tre species1line1_4971.mcc2.tre

########################################A lineage的种内分歧
*****getAline.sh
#!/bin/bash
i=0
cat singlecopygene1_1_more1000_assembly_addHv2_4971.txt | while read line
do
    i=$((i+1))
    echo $i
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$5}' > nameA.txt
		nameA=$(sed -n 1p nameA.txt)
    echo $nameA
		cat ./A_alin/subA_${nameA}.afa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > ./A_alin/A_1line_$i.alin.fa
		sed -i '1d' ./A_alin/A_1line_$i.alin.fa
done
**产生.xml文件
#!/bin/bash
for i in {1..4971}
do
  java -jar /data2/xuebo/Projects/Speciation/javaCode/C39_get_xmlfile_forbeast3_A.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A_alin/A_1line_$i.alin.fa \
  --file2 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/A_1line_4914.xml --out ./xml_A/A_1line_${i}.xml
done
#!/bin/bash
for i in {1..4971}
do
	/data1/home/xuebo/software/beast/bin/beast -threads 10 /data1/home/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/xml_A/A_1line_${i}.xml &
	monitor java 22 2s
done
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 1000000 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_species/A_${i}.trees A_${i}.mcc.tre &
	monitor java 24 2s
done
********runtreeannotator.sh
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 20 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_A/A_1line_${i}.trees A_1line_${i}.mcc.tre &
	monitor java 24 2s
done
*****改名字
#!/bin/bash
for i in {1..4971}
do
    echo "TREE"${i};
		sed -n '55,55p' A_1line_${i}.mcc.tre |  sed 's/TREE1/'"TREE${i}"'/' >>  A_1line1_4971.mcc.tre
done
*******runAPE.r
#!/usr/bin/Rscript
library(ape)
setwd("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_A")
DD = read.table("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/treeannotator_A/A_1line1_4971.mcc.tre",header = F, sep = "\t",quote = "\"")
allapetree <- matrix("NA", nrow = 4971, ncol = 1)
for(i in 1:4971){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol =  1e-3,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  write.tree(res,paste("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_A/apeTree/",i,".tree",sep=""))
}
cat *tree > ../apeA_1line1_4971.mcc.tre
shuf -n 1500 apeA_1line1_4971.mcc.tre > apeA_shuf1500_1line1_4971.mcc.tre
**算时间
/data1/home/xuebo/software/beast/bin/treeannotator -burnin 10 -heights median  A_1line1_4971.mcc.tre A_1line1_4971.mc2.tre
*****getBline.sh
#!/bin/bash
i=0
cat singlecopygene1_1_more1000_assembly_addHv2_4971.txt | while read line
do
    i=$((i+1))
    echo $i
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$6}' > nameB.txt
		nameB=$(sed -n 1p nameB.txt)
    echo $nameB
		cat ./B_alin/subB_${nameB}.afa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > ./B_alin/B_1line_$i.alin.fa
		sed -i '1d' ./B_alin/B_1line_$i.alin.fa
done
**产生.xml文件
#!/bin/bash
for i in {1..4971}
do
  java -jar /data2/xuebo/Projects/Speciation/javaCode/C39_get_xmlfile_forbeast3_B.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B_alin/B_1line_$i.alin.fa \
  --file2 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/B_1line_4970.xml --out ./xml_B/B_1line_${i}.xml
done
#!/bin/bash
for i in {1..4971}
do
	/data1/home/xuebo/software/beast/bin/beast -threads 10 /data1/home/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/xml_B/B_1line_${i}.xml &
	monitor java 22 2s
done
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 20 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_B/B_1line_${i}.trees B_1line_${i}.mcc.tre &
	monitor java 24 2s
done
*****改名字
#!/bin/bash
for i in {1..4971}
do
    echo "TREE"${i};
		sed -n '51,51p' B_1line_${i}.mcc.tre |  sed 's/TREE1/'"TREE${i}"'/' >>  B_1line1_4971.mcc.tre
done
*******runAPE.r
#!/usr/bin/Rscript
library(ape)
setwd("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_B")
DD = read.table("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/treeannotator_B/B_1line1_4971.mcc.tre",header = F, sep = "\t",quote = "\"")
allapetree <- matrix("NA", nrow = 4971, ncol = 1)
for(i in 1:4971){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol =  1e-3,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  write.tree(res,paste("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_B/apeTree/",i,".tree",sep=""))
}
cat *tree > ../apeB_1line1_4971.mcc.tre
shuf -n 1500 apeB_1line1_4971.mcc.tre > apeB_shuf1500_1line1_4971.mcc.tre
**算时间
/data1/home/xuebo/software/beast/bin/treeannotator -burnin 10 -heights median  B_1line1_4971.mcc.tre B_1line1_4971.mc2.tre
*****getDline.sh
#!/bin/bash
i=0
cat singlecopygene1_1_more1000_assembly_addHv2_4971.txt | while read line
do
    i=$((i+1))
    echo $i
    echo $line > test111.txt     
		cat test111.txt   | awk '{print$7}' > nameD.txt
		nameD=$(sed -n 1p nameD.txt)
    echo $nameD
		cat ./D_alin/subD_${nameD}.afa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  > ./D_alin/D_1line_$i.alin.fa
		sed -i '1d' ./D_alin/D_1line_$i.alin.fa
done
**产生.xml文件
#!/bin/bash
for i in {1..4971}
do
  java -jar /data2/xuebo/Projects/Speciation/javaCode/C39_get_xmlfile_forbeast3_D.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D_alin/D_1line_$i.alin.fa \
  --file2 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/D_1line_4963.xml --out ./xml_D/D_1line_${i}.xml
done
#!/bin/bash
for i in {1..4971}
do
	/data1/home/xuebo/software/beast/bin/beast -threads 10 /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/xml_D/D_1line_${i}.xml &
	monitor java 22 2s
done
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 1000000 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_species/A_${i}.trees A_${i}.mcc.tre &
	monitor java 24 2s
done
#!/bin/bash
for i in {1..4971}
do
  /data1/home/xuebo/software/beast/bin/treeannotator -burnin 20 -heights median  /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/beasttree_D/D_1line_${i}.trees D_1line_${i}.mcc.tre &
	monitor java 24 2s
done
*****改名字
#!/bin/bash
for i in {1..4971}
do
    echo "TREE"${i};
		sed -n '37,37p' D_1line_${i}.mcc.tre |  sed 's/TREE1/'"TREE${i}"'/' >>  D_1line1_4971.mcc.tre
done
*******runAPE.r
#!/usr/bin/Rscript
library(ape)
setwd("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_D")
DD = read.table("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/treeannotator_D/D_1line1_4971.mcc.tre",header = F, sep = "\t",quote = "\"")
allapetree <- matrix("NA", nrow = 4971, ncol = 1)
for(i in 1:4971){
  DD.phy <- read.tree(text = as.character(DD[i,]))
  res = chronopl(DD.phy, lambda=0, age.min = 1, age.max = NULL,
                 node = "root", S = 1, tol =  1e-3,
                 CV = FALSE, eval.max = 500, iter.max = 500)
  write.tree(res,paste("/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/beast/ape_D/apeTree/",i,".tree",sep=""))
}
cat *tree > ../apeD_1line1_4971.mcc.tre
shuf -n 1500 apeD_1line1_4971.mcc.tre > apeD_shuf1500_1line1_4971.mcc.tre
**算时间
/data1/home/xuebo/software/beast/bin/treeannotator -burnin 10 -heights median  D_1line1_4971.mcc.tre D_1line1_4971.mc2.tre

#############使用这4971个基因做树
#/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/tree_4971/nobarley
vcf-concat chr1.vcf chr2.vcf chr7.vcf chr8.vcf chr13.vcf chr14.vcf chr19.vcf chr20.vcf chr25.vcf chr26.vcf chr31.vcf chr32.vcf chr37.vcf chr38.vcf > ./A/Alineage_4971.vcf
279,240
vcf-concat chr3.vcf chr4.vcf chr9.vcf chr10.vcf chr15.vcf chr16.vcf chr21.vcf chr22.vcf chr27.vcf chr28.vcf chr33.vcf chr34.vcf chr39.vcf chr40.vcf > ./B/Blineage_4971.vcf
378,976
vcf-concat chr5.vcf chr6.vcf chr11.vcf chr12.vcf chr17.vcf chr18.vcf chr23.vcf chr24.vcf chr29.vcf chr30.vcf chr35.vcf chr36.vcf chr41.vcf chr42.vcf > ./D/Dlineage_4971.vcf
92,863
***A
bgzip -c Alineage_4971.vcf > Alineage_4971.vcf.gz
tabix Alineage_4971.vcf.gz
java -jar /data1/home/xuebo/software/WGS.jar --model vcf --type toFasta --file Alineage_4971.vcf.gz --out Alineage_4971.fasta &
java -jar /data2/xuebo/Projects/Speciation/javaCode/C23_fasta2phy.jar --file1 Alineage_4971.fasta --out Alineage_4971.phy &
nohup raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s Alineage_4971.phy -n Alineage_4971.raxml -o D014 -T 50 > nohupAlineage_4971 2>& 1 &
#/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/tree_4971/withbarley
vcf-concat chr1.withBarley.vcf chr2.withBarley.vcf chr7.withBarley.vcf chr8.withBarley.vcf chr13.withBarley.vcf chr14.withBarley.vcf chr19.withBarley.vcf chr20.withBarley.vcf chr25.withBarley.vcf chr26.withBarley.vcf chr31.withBarley.vcf chr32.withBarley.vcf chr37.withBarley.vcf chr38.withBarley.vcf > ./A/Alineage_withBarley_4971.vcf
128,269
vcf-concat chr3.withBarley.vcf chr4.withBarley.vcf chr9.withBarley.vcf chr10.withBarley.vcf chr15.withBarley.vcf chr16.withBarley.vcf chr21.withBarley.vcf chr22.withBarley.vcf chr27.withBarley.vcf chr28.withBarley.vcf chr33.withBarley.vcf chr34.withBarley.vcf chr39.withBarley.vcf chr40.withBarley.vcf > ./B/Blineage_withBarley_4971.vcf
180,486
vcf-concat chr5.withBarley.vcf chr6.withBarley.vcf chr11.withBarley.vcf chr12.withBarley.vcf chr17.withBarley.vcf chr18.withBarley.vcf chr23.withBarley.vcf chr24.withBarley.vcf chr29.withBarley.vcf chr30.withBarley.vcf chr35.withBarley.vcf chr36.withBarley.vcf chr41.withBarley.vcf chr42.withBarley.vcf > ./D/Dlineage_withBarley_4971.vcf
43,259
***A
bgzip -c Alineage_withBarley_4971.vcf > Alineage_withBarley_4971.vcf.gz
tabix Alineage_withBarley_4971.vcf.gz
java -jar /data1/home/xuebo/software/WGS.jar --model vcf --type toFasta --file Alineage_withBarley_4971.vcf.gz --out Alineage_withBarley_4971.fasta &
java -jar /data2/xuebo/Projects/Speciation/javaCode/C23_fasta2phy.jar --file1 Alineage_withBarley_4971.fasta --out Alineage_withBarley_4971.phy &
nohup raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s Alineage_withBarley_4971.phy -n Alineage_withBarley_4971.raxml -o barley -T 80 > nohupAlineage_4971 2>& 1 &
***B
bgzip -c Blineage_withBarley_4971.vcf > Blineage_withBarley_4971.vcf.gz
tabix Blineage_withBarley_4971.vcf.gz
java -jar /data1/home/xuebo/software/WGS.jar --model vcf --type toFasta --file Blineage_withBarley_4971.vcf.gz --out Blineage_withBarley_4971.fasta &
java -jar /data2/xuebo/Projects/Speciation/javaCode/C23_fasta2phy.jar --file1 Blineage_withBarley_4971.fasta --out Blineage_withBarley_4971.phy &
nohup raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s Blineage_withBarley_4971.phy -n Blineage_withBarley_4971.raxml -o barley -T 80 > nohupBlineage_4971 2>& 1 &
***D
bgzip -c Dlineage_withBarley_4971.vcf > Dlineage_withBarley_4971.vcf.gz
tabix Dlineage_withBarley_4971.vcf.gz
java -jar /data1/home/xuebo/software/WGS.jar --model vcf --type toFasta --file Dlineage_withBarley_4971.vcf.gz --out Dlineage_withBarley_4971.fasta &
java -jar /data2/xuebo/Projects/Speciation/javaCode/C23_fasta2phy.jar --file1 Dlineage_withBarley_4971.fasta --out Dlineage_withBarley_4971.phy &
nohup raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s Dlineage_withBarley_4971.phy -n Dlineage_withBarley_4971.raxml -o barley -T 80 > nohupDlineage_4971 2>& 1 &





















