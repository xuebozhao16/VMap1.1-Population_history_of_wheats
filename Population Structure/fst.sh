##算fst
vcftools --gzvcf ../ori/A.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out A.cultivar_small_landraces_small_1 &
vcftools --gzvcf ../ori/B.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out B.cultivar_small_landraces_small_1 &
vcftools --gzvcf ../ori/D.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out D.cultivar_small_landraces_small_1 &

vcftools --gzvcf ../ori/A.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_2.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out A.cultivar_small_landraces_small_2 &
vcftools --gzvcf ../ori/B.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_2.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out B.cultivar_small_landraces_small_2 &
vcftools --gzvcf ../ori/D.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_2.txt --weir-fst-pop ../group/group_cultivar_small.txt   --fst-window-size  10000 --out D.cultivar_small_landraces_small_2 &

vcftools --gzvcf ../ori/A.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_landraces_small_2.txt   --fst-window-size  10000 --out A.landraces_small_12 &
vcftools --gzvcf ../ori/B.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_landraces_small_2.txt   --fst-window-size  10000 --out B.landraces_small_12 &
vcftools --gzvcf ../ori/D.lineage.vcf.gz --weir-fst-pop ../group/group_landraces_small_1.txt --weir-fst-pop ../group/group_landraces_small_2.txt   --fst-window-size  10000 --out D.landraces_small_12 &


#############permutation test
########################################## Weinkorn & Deinkorn
#!/bin/bash
for m in {1..10}
do
  mkdir "p"$m
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/A.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_Einkorn/"p"$m".Weinkorn.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_Einkorn/"p"$m".Deinkorn.txt" --fst-window-size  1000000 --out "A_p"$m &
  cd ..
done
########################################## Wemmer & Demmer
#!/bin/bash
for m in {1..10}
do
  mkdir "p"$m
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/A.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_emmer/"p"$m".Wemmer.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_emmer/"p"$m".Demmer.txt" --fst-window-size  1000000 --out "A_p"$m &
  cd ..
done
#!/bin/bash
for m in {1..10}
do
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/B.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_emmer/"p"$m".Wemmer.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/WD_emmer/"p"$m".Demmer.txt" --fst-window-size  1000000 --out "B_p"$m &
  cd ..
done
########################################## Demmer & durum
#!/bin/bash
for m in {1..10}
do
  mkdir "p"$m
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/A.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/D_durum/"p"$m".Demmer.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/D_durum/"p"$m".durum.txt" --fst-window-size  1000000 --out "A_p"$m &
  cd ..
done
#!/bin/bash
for m in {1..10}
do
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/B.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/D_durum/"p"$m".Demmer.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/D_durum/"p"$m".durum.txt" --fst-window-size  1000000 --out "B_p"$m &
  cd ..
done
########################################## land1 & cul3
#!/bin/bash
for m in {1..10}
do
  mkdir "p"$m
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/A.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".land1.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".cul3.txt" --fst-window-size  1000000 --out "A_p"$m &
  cd ..
done
#!/bin/bash
for m in {1..10}
do
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/B.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".land1.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".cul3.txt" --fst-window-size  1000000 --out "B_p"$m &
  cd ..
done
#!/bin/bash
for m in {1..10}
do
  cd "p"$m
     vcftools --gzvcf /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/ori/D.lineage.vcf.gz --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".land1.txt" --weir-fst-pop /data1/home/yaozhou/Projects/EVO/data/lineage/final/V8/permutationFst/group/land1_cul3/"p"$m".cul3.txt" --fst-window-size  1000000 --out "D_p"$m &
  cd ..
done




################################################################permutation test 数据整理
###########WD_Einkorn
#!/bin/bash
for m in {1..10}
do
  cd "p"$m 
     sed -i '1d' "A_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "A_p"$m".windowed.weir.fst" --pos 4  --out "p"$m"_max.txt" 
  cp "p"$m"_max.txt"  ../
  cd ..
done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
###########WD_emmer
#!/bin/bash
for m in {1..10}
do
  cd "p"$m 
     sed -i '1d' "A_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "A_p"$m".windowed.weir.fst" --pos 4  --out "Ap"$m"_max.txt" 
     sed -i '1d' "B_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "B_p"$m".windowed.weir.fst" --pos 4  --out "Bp"$m"_max.txt" 
     cat "Ap"$m"_max.txt" "Bp"$m"_max.txt" > "p"$m"_max.txt" 
  cp "p"$m"_max.txt"  ../
  cd ..
done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
##merge数据			 
#!/bin/bash
		for m in {1..10}
				do
						cd "p"$m 
						rm "Ap"$m"_max.txt" 
						rm "Bp"$m"_max.txt" 
						rm "p"$m"_max.txt" 
					  cat "A_p"$m".windowed.weir.fst"  "B_p"$m".windowed.weir.fst" > "WD_emmerp"$m".fst"
					  WGS --model file --type getMax --file "WD_emmerp"$m".fst" --pos 4  --out "p"$m"_max.txt" 
						cp "p"$m"_max.txt"   ../
						cd ..
				done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
###########D_durum
#!/bin/bash
for m in {1..10}
do
  cd "p"$m 
     sed -i '1d' "A_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "A_p"$m".windowed.weir.fst" --pos 4  --out "Ap"$m"_max.txt" 
     sed -i '1d' "B_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "B_p"$m".windowed.weir.fst" --pos 4  --out "Bp"$m"_max.txt" 
     cat "Ap"$m"_max.txt" "Bp"$m"_max.txt" > "p"$m"_max.txt" 
  cp "p"$m"_max.txt"  ../
  cd ..
done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
##merge数据			 
#!/bin/bash
		for m in {1..10}
				do
						cd "p"$m 
						rm "Ap"$m"_max.txt" 
						rm "Bp"$m"_max.txt" 
						rm "p"$m"_max.txt" 
					  cat "A_p"$m".windowed.weir.fst"  "B_p"$m".windowed.weir.fst" > "D_dup"$m".fst"
					  WGS --model file --type getMax --file "D_dup"$m".fst" --pos 4  --out "p"$m"_max.txt" 
						cp "p"$m"_max.txt"   ../
						cd ..
				done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
#############land1_cul3
#!/bin/bash
for m in {1..10}
do
  cd "p"$m 
     sed -i '1d' "A_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "A_p"$m".windowed.weir.fst" --pos 4  --out "Ap"$m"_max.txt" 
     sed -i '1d' "B_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "B_p"$m".windowed.weir.fst" --pos 4  --out "Bp"$m"_max.txt" 
     sed -i '1d' "D_p"$m".windowed.weir.fst" 
     WGS --model file --type getMax --file "D_p"$m".windowed.weir.fst" --pos 4  --out "Dp"$m"_max.txt" 
     cat "Ap"$m"_max.txt" "Bp"$m"_max.txt" "Dp"$m"_max.txt"  > "p"$m"_max.txt" 
  cp "p"$m"_max.txt"  ../
  cd ..
done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt
##merge数据			 
#!/bin/bash
		for m in {1..10}
				do
						cd "p"$m 
						rm "Ap"$m"_max.txt" 
						rm "Bp"$m"_max.txt" 
						rm "Dp"$m"_max.txt"
						rm "p"$m"_max.txt" 
					  cat "A_p"$m".windowed.weir.fst" "B_p"$m".windowed.weir.fst" "D_p"$m".windowed.weir.fst" > "land1_cul3p"$m".fst"
					  WGS --model file --type getMax --file "land1_cul3p"$m".fst" --pos 4  --out "p"$m"_max.txt" 
						cp "p"$m"_max.txt"   ../
						cd ..
				done
cat p*_max.txt > pAll_max.txt
WGS --model file --type getMax --file pAll_max.txt --pos 2 --out A_centi.txt
WGS --model file --type getMax --file pAll_max.txt --pos 3 --out A_milli.txt


###############sort & 按照阈值取相应的XPCLR的值
## 加-r从大到小排序 -k是对哪一列进行排序 -g是对科学计数法排序
cp wild_dom_einkorn.10K.windowed.weir.fst ./sortfst/A.fst.txt
cat A.wild_dom_emmer.10K.windowed.weir.fst B.wild_dom_emmer.10K.windowed.weir.fst > ./sortfst/AB_15.fst.txt
cat A.dom_durum.10K.windowed.weir.fst B.dom_durum.10K.windowed.weir.fst > ./sortfst/AB_45.fst.txt
cat A.cultivar_small_landraces_small_1.windowed.weir.fst B.cultivar_small_landraces_small_1.windowed.weir.fst D.cultivar_small_landraces_small_1.windowed.weir.fst > ./sortfst/ABD.fst.txt

sort -k6 -rg A.fst.txt > A.Sortfst.txt
sort -k6 -rg AB_15.fst.txt  > AB_15.Sortfst.txt
sort -k5 -rg AB_45.fst.txt  > AB_45.Sortfst.txt
sort -k5 -rg ABD.fst.txt > ABD.Sortfst.txt

/data1/home/yaozhou/Projects/EVO/data/lineage/final/V10/Fst/10K
###############sort & 按照阈值取相应的XPCLR的值
cp lineageA.1M.WDeinkorn.bed.fst A.fst.txt
cat lineageA.1M.WD.emmer.bed.fst lineageB.1M.WDemmer.bed.fst > AB_15.fst.txt
cat lineageA.1M.DD.bed.fst lineageB.1M.DD.bed.fst > AB_45.fst.txt
cat A.cultivar_small_landraces_small_1.windowed.weir.fst B.cultivar_small_landraces_small_1.windowed.weir.fst D.cultivar_small_landraces_small_1.windowed.weir.fst > ./sortfst/ABD.fst.txt

ForManhattanPlot.java

cat Chr_A.fst.txt |grep chr1A > fst_chr1A_Chr_A.txt
cat Chr_AB_15.fst.txt |grep chr1A > fst_chr1A_Chr_AB_15.txt
cat Chr_AB_45.fst.txt |grep chr1A > fst_chr1A_Chr_AB_45.txt
cat Chr_AB_cul5.fst.txt |grep chr1A > fst_chr1A_Chr_AB_cul5.txt



###计算reference allele frequency
/data1/home/yaozhou/Projects/EVO/data/lineage/final/V11/XPCLR_group
java -jar -Xms10g -Xmx30g ../Wheat.jar --file ../Alineage.txt  --out Alineage_1000snp.txt &
#!/bin/bash
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
     java -jar -Xms10g -Xmx30g ../../Wheat.jar  --file "Agroup4Chr"$i".vcf.gz" --out "Allele_freq_Agroup4Chr"$i".txt"  &
     COUNT=$(ps -ef |grep "java" |wc -l)
			 echo $COUNT
			 while [ $COUNT -gt 25 ]
			 do
		       COUNT=$(ps -ef |grep "java" |wc -l)
		       sleep 60s
		   done
done

#!/bin/bash
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
     java -jar -Xms10g -Xmx30g ../../Wheat.jar  --file "ABgroupculChr"$i".vcf.gz" --out "Allele_freq_ABgroupculChr"$i".txt"  &
     COUNT=$(ps -ef |grep "java" |wc -l)
			 echo $COUNT
			 while [ $COUNT -gt 25 ]
			 do
		       COUNT=$(ps -ef |grep "java" |wc -l)
		       sleep 60s
		   done
done
#!/bin/bash
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
    WGS --chr $i --model diversity --type bedPi --file "Allele_freq_ABgroup1Chr"$i".txt" --file2 /data1/home/yaozhou/Projects/EVO/data/merge/lineage/bed/all.bed2  --out "ABgroup1_10k_Chr"$i".txt" &
done
for ((i=1;i<=42;i++))do echo ABgroup1_10k_Chr$i.txt;done | xargs -i cat {} >> Allele_freq_ABgroup1_10k.txt




###fst和maf在一起画图
这是去掉nan
grep -v "nan" WDEI.fst > WDEI.clean.fst
grep -v "nan" WDEM.fst > WDEM.clean.fst
grep -v "nan" DD.fst > DD.clean.fst
grep -v "nan" LAN_CUL.fst > LAN_CUL.clean.fst


































