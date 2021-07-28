######现在是对hanplotype block进行分析
******先做个小的测试
bcftools view -r 5:61000000-77000000 /data2/xuebo/Projects/Speciation/E3/chr5.all.vcf.gz  -o Haplo_chr5_61_77M.vcf
bcftools view -r 5:104000000-105000000 /data2/xuebo/Projects/Speciation/E3/chr5.all.vcf.gz  -o Haplo_chr5_104_105M.vcf
bcftools view  -S /data2/xuebo/Projects/Speciation/group/AABBDD_taxa_NoSynthetic.txt Haplo_chr5_104_105M.vcf --threads 4 -o Haplo_chr5_104_105M_ABD.vcf
**本地java
String infileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block/Haplo_chr5_104_105M_ABD.vcf";
String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/lostruct/test_block/Haplo_chr5_104_105M_ABD.txt";
new ForHeatmap(infileS,outfileS);











