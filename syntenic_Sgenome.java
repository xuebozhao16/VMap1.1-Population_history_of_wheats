/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import static EvolutionWheat.ForVcftoolsGroup.getTextReader;
import static EvolutionWheat.ForVcftoolsGroup.listFilesEndsWith;
import static EvolutionWheat.ForVcftoolsGroup.listRecursiveFiles;
import com.google.common.collect.Sets;
import com.google.common.collect.Sets.SetView;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.IntStream;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
//import static utils.IOUtils.getTextReader;
//import static utils.IOUtils.listFilesEndsWith;
//import static utils.IOUtils.listRecursiveFiles;
//import utils.PArrayUtils;

/**
 *
 * @author xuebozhao
 */
public class syntenic_Sgenome {
//    public syntenic_Sgenome(String infileS1,String infileS2,String outfileS){
//        //this.getSgenome_SyntenicSite(infileS1, infileS2, outfileS);
//        //this.getSgenome_SyntenicVCF(infileS1, infileS2, outfileS);
//       
//        //this.getTwoVCFfile_merge(infileS1, infileS2, outfileS);
//        
//    }
    public syntenic_Sgenome(){
        //this.mergePosList();
    }
//    public syntenic_Sgenome(String infileS,String outfileS,String ref){
//        //this.readfilder(infileS, outfileS, ref);
//    }
    
    public syntenic_Sgenome(String infileS1,String outfileS1,String outfileS2){
        //this.getScan_hapPosANDposAllelefiles(infileS1, outfileS1, outfileS2);
        this.C19_getScan_hapPosANDposAllelefiles(infileS1, outfileS1, outfileS2);
          
    }
    
    public syntenic_Sgenome(String infileS,String outfileS){
        this.getNewTaxafile(infileS, outfileS);
        //this.getVCFfile_info1(infileS, outfileS);  
        //this.C3_DepthandSD(infileS, outfileS);
        //this.C4_DepthandSD(infileS, outfileS);
        //this.getVCFfile_info4(infileS, outfileS);
    }
   
    public syntenic_Sgenome(String infileS1,String infileS2,String infileS3,String outfileS){
        //this.getVCFfile_info2(infileS1, infileS2, infileS3, outfileS);
        this.getVCFfile_info3(infileS1, infileS2, infileS3, outfileS);
    }
    
    
    //这个方法是输入depth文件和B基因组所有的syntenic site的文件，得到的是depth在40-150(40-120)之间的而且是syntenic的site
    //java -jar -Xms10g -Xmx50g /data1/home/xuebo/Projects/Speciation/javaCode/getSgenome_SyntenicSite.jar --file1 /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/all_SyntenicSite/chr${chr}.txt \
  //--file2 /data1/home/xuebo/Projects/Speciation/bamS10/Depth/depth.${chr}.txt.gz --out /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/getSgenome_SyntenicSite/S001.chr${chr}.SyntenicSite.txt & 

    public void C2_getSgenome_SyntenicSite(String infileS1,String infileS2,String outfileS){
        try {
            String temp1 = null;
            String temp2 = null;
            BufferedReader br1;
            BufferedReader br2;
            if (infileS1.endsWith("gz")) {
                br1 = IOUtils.getTextGzipReader(infileS1);
            } else {
                br1 = IOUtils.getTextReader(infileS1);
            }
            if (infileS2.endsWith("gz")) {
                br2 = IOUtils.getTextGzipReader(infileS2);
            } else {
                br2 = IOUtils.getTextReader(infileS2);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //把基因组所有的syntenic site的文件放进一个set里面
            Set SyntenicSite = new HashSet();
            while((temp1 = br1.readLine())!= null){
                SyntenicSite.add(temp1.split("\t")[1]);
            }
            System.out.print("readed site file" + "\n");
            while((temp2 = br2.readLine())!= null){
                String tem[] = temp2.split("\t");
                int depthSum = 0;
                if(!SyntenicSite.add(tem[1])){
                    for(int i=2;i<12;i++){
                        depthSum = depthSum + Integer.valueOf(tem[i]);
                    }
                    if(depthSum<=120 && depthSum>=40){
                        bw.write(tem[0] + "\t" + tem[1] + "\t" + depthSum + "\n");
                    }
                }else{
                    SyntenicSite.remove(tem[1]);
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //得到depth和SD的关系
    public void C3_DepthandSD(String infileS,String outfileS){
         try{    
            String temp = null; 
            //String temporder = null;
            BufferedReader Depth = null;
            //BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            HashMap<String, Double> hashMapMean = new HashMap<String, Double>();
            File f = new File(infileS);
            File[] fs = IOUtils.listRecursiveFiles(f);
            File[] sub = IOUtils.listFilesEndsWith(fs, "gz");
            for(File fi:sub){
                Depth = IOUtils.getTextGzipReader(fi.toString());
                //Depth = IOUtils.getTextReader(fi.toString());
                String taxaNamelist[] = fi.toString().split("/");
                String taxaName = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
                System.out.print("It's " + taxaName + "\n");
                while((temp = Depth.readLine()) != null){
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    String[] tem = temp.split("\t");
                    //int sum = 0;
                    //double sd = 0;             
                    for(int i=2;i<12;i++){
                        //System.out.println(tem[i]+"\n");
                        depthList.add(Double.valueOf(tem[i]));
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    //System.out.println(d.getMean()+"\n");
                    bw.write(d.getMean() + "\t");
                    bw.write(d.getStandardDeviation() + "\n"); 
                }
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    //得到depth和SD的关系,这次是分着染色体的
        public void C4_DepthandSD(String infileS,String outfileS){
             try{    
                String temp = null; 
                //String temporder = null;
                BufferedReader br = null;
                if (infileS.endsWith("gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    } else {
                        br = IOUtils.getTextReader(infileS);
                    }
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                while((temp = br.readLine()) != null){
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    String[] tem = temp.split("\t");
                    //int sum = 0;
                    //double sd = 0;             
                    for(int i=2;i<12;i++){
                        //System.out.println(tem[i]+"\n");
                        depthList.add(Double.valueOf(tem[i]));
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    //System.out.println(d.getMean()+"\n");
                    bw.write(d.getMean() + "\t");
                    bw.write(d.getStandardDeviation() + "\n"); 
                }        
                bw.flush();
                bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
            }
        }
    
    //这个方法是为了得到符合depth分布和syntenic分布的VCF文件
    // java -jar -Xms40g -Xmx70g /data1/home/xuebo/Projects/Speciation/javaCode/getSvcf_SyntenicSite.jar --file1 /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/getSgenome_SyntenicSite/S.chr${chr}.SyntenicSite_120.txt \
    //--file2 /data1/home/xuebo/Projects/Speciation/gatkS10/combinegvcf/Sgenome.chr${chr}.vcf.gz --out /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/vcf_SyntenicSite/S.chr${chr}.vcf.gz & 

    public void getSgenome_SyntenicVCF(String infileS1,String infileS2,String outfileS){
        try {
            String temp1 = null;
            String temp2 = null;
            BufferedReader br1;
            BufferedReader br2;
            if (infileS1.endsWith("gz")) {
                br1 = IOUtils.getTextGzipReader(infileS1);
            } else {
                br1 = IOUtils.getTextReader(infileS1);
            }
            if (infileS2.endsWith("gz")) {
                br2 = IOUtils.getTextGzipReader(infileS2);
            } else {
                br2 = IOUtils.getTextReader(infileS2);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            //把基因组所有的syntenic site的文件放进一个set里面
            Set SyntenicSite = new HashSet();
            while((temp1 = br1.readLine())!= null){
                SyntenicSite.add(temp1.split("\t")[1]);
            }
            System.out.print("readed SyntenicSite file" + "\n");
            //现在开始读VCF文件
            while((temp2 = br2.readLine())!= null){
                if(temp2.startsWith("#")){
                    bw.write(temp2 + "\n");
                }else{
                   String tem[] = temp2.split("\t"); 
                   if(!SyntenicSite.add(tem[1])){
                       bw.write(temp2 + "\n");
                   }else{
                    SyntenicSite.remove(tem[1]);
                   }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    
    //这个方法是为了把两个vcf文件合并，S和V11的Blineage合成一个vcf文件，使用的是V11的表头
    //	java -jar -Xms100g -Xmx120g /data1/home/xuebo/Projects/Speciation/javaCode/getTwoVCFfile_merge.jar --file1 /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/vcf_SyntenicSite/S_biallelic.chr${chr}.recode.vcf.gz \
	//--file2 /data1/home/xuebo/Projects/Evo/xpclr/lineage/V11/chr${chr}.vcf.gz --out /data1/home/xuebo/Projects/Speciation/syntenic_Sgenome/mergedVCF_withV11/E1.chr${chr}.vcf & 

    public void getTwoVCFfile_merge(String infileS1,String infileS2,String outfileS){
        try{
            String temp = null;
            String temp2 = null;
            int Alength = 0;
            int Blength = 0;
            int sum = 0;
            StringBuilder headlineSB1 = new StringBuilder();
            StringBuilder headlineSB2 = new StringBuilder();
            HashMap<Integer, String> hashMap1 = new HashMap<Integer, String>();
            BufferedReader br1;
            BufferedReader br2;
            if (infileS1.endsWith("gz")) {
                br1 = IOUtils.getTextGzipReader(infileS1);
            } else {
                br1 = IOUtils.getTextReader(infileS1);
            }
            if (infileS2.endsWith("gz")) {
                br2 = IOUtils.getTextGzipReader(infileS2);
            } else {
                br2 = IOUtils.getTextReader(infileS2);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //现在读的是S genome的VCF文件
            while((temp = br1.readLine()) != null){
                String tem1[] = temp.split("\t");
                Alength = tem1.length - 9;
                //String prefix = tem1[0].substring(0, 2);
                if(tem1[0].startsWith("##")){
                    //前边的不要
                }
                else if(tem1[0].startsWith("#C")){
                    //String headline = br.readLine();
                    String headtem[] = temp.split("\t");
                    for(int j=0; j<headtem.length;j++){
                        //这里是保存一下S genome 的vcf的表头 S001   S002.....S010
                        headlineSB1.append(headtem[j]).append("\t");
                    }
                }
                else {
                    hashMap1.put(Integer.valueOf(tem1[1]), temp);
                }
            }
            System.out.println("readed the Sgenome VCF file" + "\n");
            //现在开始读V11的B lineage
            while((temp2 = br2.readLine()) != null){
                String tem[] = temp2.split("\t");
                Blength = tem.length - 9;
                //String prefix = tem[0].substring(0, 2);
                if(tem[0].startsWith("##")){
                    bw.write(temp2 + "\n");  //因为想要V11的表头，所以这个前边的是输出的
                }
                else if(tem[0].startsWith("#C")){
                    //String headline = br.readLine();
                    String headtem[] = temp2.split("\t");
                    for(int j=9; j<headtem.length;j++){
                        //这里是保存V11 Blineage的vcf的表头
                        headlineSB2.append(headtem[j]).append("\t");
                    }
                    //System.out.println(headlineSB1.toString() + headlineSB2.toString() + "\n");
                    bw.write(headlineSB1.toString() + headlineSB2.toString() + "\n");  //现在是把两个表头一起输出来
                }
                else {
                    int PosAB = Integer.valueOf(tem[1]);
                    if(hashMap1.get(PosAB) == null){  //这里指的是S genome里面没有V11 Blineage的这个位点,V11 Blineage的这个点是独有的
                        StringBuilder headlineABhead = new StringBuilder();
                        StringBuilder headlineAB = new StringBuilder();
                        for(int i=0;i<8;i++){
                            headlineABhead.append(tem[i] + "\t");
                        }
                        headlineABhead.append(tem[8]);  //把这个位点的V11 Blineage的前八列放到headlineABhead里面
                        for(int i=9;i<(tem.length-1);i++){
                            headlineAB.append(tem[i] + "\t");
                        }
                        headlineAB.append(tem[tem.length-1]); //把这个位点的V11 Blineage的第九列之后的放到headlineAB里面
                        String valueNew = headlineABhead.toString() + "\t" + get00(Alength) + "\t" + headlineAB.toString(); //把S genome的位置用0/0填充起来
                        hashMap1.put(Integer.valueOf(tem[1]), valueNew);
                    }else{  //这里指的是S genome里面有V11 Blineage的这个位点
                        StringBuilder headlineAB = new StringBuilder(); 
                        //现在是判断两个ALT是否一样，要是不一样的话就去掉
                        String ABalt = tem[4];
                        //System.out.println(tem[1] + "\t" + ABalt);
                        String Aalt = hashMap1.get(Integer.valueOf(tem[1])).split("\t")[4];
                        //System.out.println(Aalt);
                        if(ABalt.equals(Aalt)){
                            for(int i=9;i<(tem.length-1);i++){
                            headlineAB.append(tem[i] + "\t");
                            }
                            headlineAB.append(tem[tem.length-1]);  //把这个位点的V11 Blineage的第九列之后的放到headlineAB里面
                            String valueNew = hashMap1.get(Integer.valueOf(tem[1])) + "\t" + headlineAB.toString();  //S genome里面有V11 Blineage的这个位点,所以两个部分直接合起来
                            hashMap1.replace(Integer.valueOf(tem[1]), valueNew);
                        }else{
                            hashMap1.remove(Integer.valueOf(tem[1]), hashMap1.get(Integer.valueOf(tem[1])));
                            sum = sum +1;
                        }
             }                   
                }
            }
            System.out.println(sum);
            //hashMap按照key值进行排序
            Object[] keyall = hashMap1.keySet().toArray();   
            Arrays.sort(keyall);  
            for (int n = 0; n < keyall.length; n++) {   
                String valueall[] = hashMap1.get(keyall[n]).split("\t");
                if(valueall.length == Alength + Blength + 9){ //现在是指达到了那个长度，即S+Blineage+表头的长度
                    bw.write(hashMap1.get(keyall[n]) + "\n");
                }else {
                    String valueNew = hashMap1.get(keyall[n]) + "\t" + get00(Blength); //没有S+Blineage+表头的长度，表示只有S基因组有，但是Blineage没有，需要用0/0填满
                    //System.out.println(valueNew);
                    bw.write(valueNew + "\n");
                }               
            }   
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    // 这个方法是产生指定个数的NA
    private String getNA(Integer Alength){
        StringBuilder StringNA = new StringBuilder();
        for(int j = 0; j< Alength-1; j++){
            StringNA.append("NA" + "\t");
        }
        StringNA.append("NA");
        return StringNA.toString();
    }
    
    
    // 这个方法是产生指定个数的0/0
    private String get00(Integer Alength){
        StringBuilder String00 = new StringBuilder();
        for(int j = 0; j< Alength-1; j++){
            String00.append("0/0" + "\t");
        }
        String00.append("0/0");
        return String00.toString();
    }
    
    
    
    //这是鲁老师写的合并两个vcf的代码,这样是得出所有位点的"Chr\tPos\tRef\tAlt"
    public void mergePosList () {
        String inFileS1 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/vcfE1/test/testsite11.vcf.gz";
        String inFileS2 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/vcfE1/test/testBlineag.vcf.gz";
        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/vcfE1/mergelist.txt";
        String[] alleles = {"A", "C", "G", "T", "D", "I"};
        Arrays.sort(alleles);
        double[] fre = new double[6];
        try {
            int chr1 = Integer.MIN_VALUE;
            int chr2 = Integer.MIN_VALUE;
            int taxaNum1 = Integer.MIN_VALUE;
            int taxaNum2 = Integer.MIN_VALUE;
            TIntArrayList posList1 = new TIntArrayList();
            List<String> referList1 = new ArrayList<>();
            List<String> altList1 = new ArrayList<>();
            List<String> altDepthList1 = new ArrayList<>();
            TIntArrayList posList2 = new TIntArrayList();
            List<String> referList2 = new ArrayList<>();
            List<String> altList2 = new ArrayList<>();
            List<String> altDepthList2 = new ArrayList<>();
            BufferedReader br = IOUtils.getTextGzipReader(inFileS1);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum1 = temp.split("\t").length-9;
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                //temp = temp.substring(0, 1000);
                tem = temp.split("\t");
                chr1 = Integer.parseInt(tem[0]);
                posList1.add(Integer.parseInt(tem[1]));
                referList1.add(tem[3]);
                altList1.add(tem[4]);
                //altDepthList1.add(tem[7].split(";")[1].replace("AD=", ""));
                altDepthList1.add(tem[7].split(";")[1].split("=")[1]);
            }
            br.close();
            br = IOUtils.getTextGzipReader(inFileS2);
            temp = null;
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum2 = temp.split("\t").length-9;
            double weight1 = (double)taxaNum1/(taxaNum1+taxaNum2);
            double weight2 = (double)taxaNum2/(taxaNum1+taxaNum2);
            tem = null;
            while ((temp = br.readLine()) != null) {
                //temp = temp.substring(0, 50);
                tem = temp.split("\t");
                chr2 = Integer.parseInt(tem[0]);
                posList2.add(Integer.parseInt(tem[1]));
                referList2.add(tem[3]);
                altList2.add(tem[4]);
                String NNN  = "0.188";
                altDepthList2.add(NNN);
            }
            if (chr1 != chr2) {
                System.out.println("Wrong input files! Program quits.");
                System.exit(0);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt");
            bw.newLine();
            TIntHashSet mergedPosSet = new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos = mergedPosSet.toArray();
            Arrays.sort(mergedPos);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < mergedPos.length; i++) {
                sb = new StringBuilder();
                int index1 = posList1.binarySearch(mergedPos[i]);
                int index2 = posList2.binarySearch(mergedPos[i]);
                if (index1 < 0 && index2 > -1) { //这个位置在posList1里面没有，在posList2里面有
                    if (altList2.get(index2).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList2.get(index2)).append("\t").append(referList2.get(index2)).append("\t").append(altList2.get(index2));
                }
                else if (index1 > -1 && index2 < 0) {  //这个位置在posList2里面没有，在posList1里面有
                    if (altList1.get(index1).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t").append(altList1.get(index1));
                }
                else  { //这个位置在posList2里面有，在posList1里面有
                    for (int j = 0; j < fre.length; j++) {
                        fre[j] = -1;
                    }
                    tem = altList1.get(index1).split(",");
                    String[] fretem = altDepthList1.get(index1).split(",");
                    double[] depth = new double[fretem.length];
                    double[] fre1 = new double[depth.length];
                    double sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        System.out.println(fretem[j]);
                        depth[j] = Double.parseDouble(fretem[j]);
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre1[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        System.out.println(fre1[j]);
                        fre[index] = fre1[j]*weight1;
                    }
                    
                    tem = altList2.get(index2).split(",");
                    fretem = altDepthList2.get(index2).split(",");
                    depth = new double[fretem.length];
                    double[] fre2 = new double[depth.length];
                    sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Double.parseDouble(fretem[j]);
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre2[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        if (fre[index] < 0) {
                            fre[index] = fre2[j+1]*weight2;
                        }
                        else {
                            fre[index] = fre[index]+fre2[j]*weight2;
                        }
                    }
                    
                    //int[] indices = PArrayUtils.getIndexByDescendingValue(fre);
                    int[] indices = PArrayUtils.getIndicesByAscendingValue(fre);
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t");
                    for (int j = 0; j < 2; j++) {
                        if (fre[indices[j]] > 0) {
                            sb.append(alleles[indices[j]]).append(",");
                        }
                        else {
                            break;
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    //现在是整理hapscan的文件 Taxa    Reference       BamPath
    //读文件夹里面的所有的以我想要的后缀(.bam)的文件
    //java -jar /data2/xuebo/Projects/Speciation/javaCode/getFileTaxa.jar  --file /data3/wgs/bam/A --out /data2/xuebo/Projects/Speciation/hapScan/TaxaRefBam/aoyue_A.txt --p /data1/home/xuebo/Projects/Speciation/reference/a_iwgscV1_forGATK.fa.gz
    public void readfilder(String infileS,String outfileS,String ref){
         try{    
            String temp = null; 
            BufferedReader xpclrFile = null;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            File f = new File(infileS);
            File[] fs = IOUtils.listRecursiveFiles(f);
            File[] sub = IOUtils.listFilesEndsWith(fs, ".rmdup.bam");
            for(File fi:sub){
                //xpclrFile = getTextReader(fi.toString());
                String taxaNamelist[] = fi.toString().split("/");
                String taxaName = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
                //String chrname = taxaName.split("_")[1];
                //String chr = chrname.split("r")[1];
                //String Blineage = "/data1/home/xuebo/Project/reference/b_iwgscV1_forGATK.fa.gz";
                bw.write(taxaName + "\t" + ref + "\t" + fi + "\n");
            } 
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //产生hapPos文件和posAllele文件
    //	java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/getScan_hapPosANDposAllelefiles.jar --file /data2/xuebo/Projects/Speciation/E1/E1.chr${chr}.vcf.gz --out /data2/xuebo/Projects/Speciation/hapScan/hapPos/chr${chr}.pos.txt \
	//--out2 /data2/xuebo/Projects/Speciation/hapScan/posAllele/chr${chr}.allele.txt &
    
    public void C19_getScan_hapPosANDposAllelefiles(String infileS1,String outfileS1,String outfileS2){
        try {
            String temp = null;
            BufferedReader br;
            if (infileS1.endsWith("gz")) {
                br = IOUtils.getTextGzipReader(infileS1);
            } else {
                br = IOUtils.getTextReader(infileS1);
            }
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            bw2.write("Chr\tPos\tRef\tAlt");
            bw2.newLine();
            //现在开始读VCF文件
            while((temp = br.readLine())!= null){
                if(temp.startsWith("#")){
                    
                }else{
                   //System.out.println(temp);
                   String tem[] = temp.split("\t");
                   bw1.write(tem[0] + "\t" + tem[1] + "\n");
                   bw2.write(tem[0] + "\t" + tem[1] + "\t" + tem[3] + "\t" + tem[4] + "\n");                  
                }
            }
            bw1.flush();
            bw2.flush();
            bw1.close();
            bw2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    
    
    
    //产生parameters文件
//        String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatScan/hapScanner2/Parameters_hapScanner2_chr1.txt";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/hapscan/parameters_file/p_DD/";
//        new syntenic_Sgenome(infileS,outfileS);
    public void getNewTaxafile(String infileS,String outfileS){
        try {            
            BufferedWriter bw = null;
            String temp = "";
            String[] te = null;
            
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38};
            //int chrA[] = {3,4,9,10,15,16,21,22,27,28,33,34,39,40};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40,5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            int chrA[] = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            //int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            
            for(int i = 0; i < chrA.length;i++ ){
                BufferedReader br = IOUtils.getTextReader(infileS);
                int line = 1;
                //bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerAABBDD_chr"+chrA[i]+".txt");
                //bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerBarley_chr"+chrA[i]+".txt");
                //bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerAe_chr"+chrA[i]+".txt");
                bw = IOUtils.getTextWriter(outfileS+"/parameters_hapScannerTu_chr"+chrA[i]+".txt");
                while((temp = br.readLine())!=null){
                    if(line==9){
                        //bw.write("/data2/xuebo/Projects/Speciation/introgression/scanAe/TaxaRefBam/TaxaRefBamChr" + chrA[i]);
                        bw.write("/data2/xuebo/Projects/Speciation/introgression/scanTu/TaxaRefBam/TaxaRefBamChr" + chrA[i]);
                        //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/TaxaRefBam/TaxaRefBamAABBDD.txt");
                        bw.newLine();
                    }else if(line==11){
                        //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/posAllele/chr"+chrA[i]+".allele.txt");
                        bw.write("/data2/xuebo/Projects/Speciation/hapScan/posAllele/chr"+chrA[i]+".allele.txt");
                        bw.newLine();
                    }
                    else if(line==13){
                        //bw.write("/data2/xuebo/Projects/Speciation/E5/hapSacn/hapPos/chr"+chrA[i]+".pos.txt");
                        bw.write("/data2/xuebo/Projects/Speciation/hapScan/hapPos/chr"+chrA[i]+".pos.txt");
                        bw.newLine();
                    }
                    else if(line==15){
                        bw.write(Integer.toString(chrA[i]));
                        bw.newLine();
                    }else if(line==17){
                        bw.write("/data1/home/xuebo/anaconda3/bin/samtools");
                        bw.newLine();
                    }else if(line==21){
                        //bw.write("/data2/xuebo/Projects/Speciation/introgression/scanAe/outchr"+chrA[i]);
                        bw.write("/data2/xuebo/Projects/Speciation/introgression/scanTu/outchr"+chrA[i]);
                        bw.newLine();
                    }else{
                        bw.write(temp);
                        bw.newLine();
                    }
                    line++;
                }
                bw.flush();
                bw.close();
            }
        }   
        catch(Exception e){
            e.printStackTrace();
        }
    }

    //现在对scan之后的E2进行SNP数量进行统计,对文件夹内的所有文件进行操作
    //这个最后没用，用的是wc -l
    public void getVCFfile_info1(String infileS,String outfileS){
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //现在开始读VCF文件
            File f = new File(infileS);
            File[] fs = IOUtils.listRecursiveFiles(f);
            File[] sub = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
            int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(sub.length, 120); //设置线程数
            for (int i = 0; i < indices.length; i++) {
                Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
                for (int j = 0; j < subLibIndices.length; j++) {
                    subLibIndices[j] = indices[i][0]+j;
                }
                List<Integer> integerList=Arrays.asList(subLibIndices);
                integerList.parallelStream()
                        .forEach((Integer index)-> {
                            try{
                                BufferedReader br = IOUtils.getTextGzipReader(sub[index].toString()); //先要读文件
                                String name[] = sub[index].toString().split("/");
                                String filename = name[name.length-1];
                                int count = 0;
                                String temp = null;
                                while((temp = br.readLine()) != null){
                                    if(temp.startsWith("#")){

                                    }else{
                                        count = count + 1;
                                    }                  
                                }
                                System.out.println(filename + "\n");
                                bw.write(filename + "\t" + count + "\n");   
                            } catch (Exception e) {
                                e.printStackTrace();
                            }   
                        });
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //在各个群体里面的分离情况进行统计
    //对Alineage来说，输入文件有3个，AA的个体的VCF文件，AABB的个体的VCF文件，AABBDD的个体的VCF文件
//    java -jar -Xms50g -Xmx70g /data2/xuebo/Projects/Speciation/javaCode/getVCFfile_info2.jar --file1 /data2/xuebo/Projects/Speciation/E2_all/macSS/chr${chr}.macSS.vcf.gz \
//	--file2 /data2/xuebo/Projects/Speciation/E2_all/macAABB/chr${chr}.macAABB.vcf.gz  --file3 /data2/xuebo/Projects/Speciation/E2_all/macAABBDD/chr${chr}.macAABBDD.vcf.gz \
//	--out /data2/xuebo/Projects/Speciation/E2_all/macSS/info2_chr${chr}.txt &
    
    public void getVCFfile_info2(String infileS1,String infileS2,String infileS3,String outfileS){
        try{
            //String temp = null;
            String temp1 = null;
            String temp2 = null;
            String temp3 = null;
            //int countall = 0;
            Set AAset = new HashSet();
            Set AABBset = new HashSet();
            Set AABBDDset = new HashSet();
            //HashMap<Integer, String> hashMap1 = new HashMap<Integer, String>();
            BufferedReader br1; //这里读的是AA个体的文件夹
            if (infileS1.endsWith("gz")) {
                br1 = IOUtils.getTextGzipReader(infileS1);
            } else {
                br1 = IOUtils.getTextReader(infileS1);
            }
            BufferedReader br2; //这里读的是AABB个体的文件夹
            if (infileS2.endsWith("gz")) {
                br2 = IOUtils.getTextGzipReader(infileS2);
            } else {
                br2 = IOUtils.getTextReader(infileS2);
            }
            BufferedReader br3; //这里读的是AABBDD个体的文件夹
            if (infileS3.endsWith("gz")) {
                br3 = IOUtils.getTextGzipReader(infileS3);
            } else {
                br3 = IOUtils.getTextReader(infileS3);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS); //写出文件，文件格式是txt格式
            //开始读genome的VCF文件
            while((temp1 = br1.readLine()) != null){
                String tem1[] = temp1.split("\t");
                if(tem1[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AAset.add(Integer.valueOf(tem1[1]));                  
                }
            }
            System.out.println("Readed AA" + "\n");
            while((temp2 = br2.readLine()) != null){
                String tem2[] = temp2.split("\t");
                if(tem2[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AABBset.add(Integer.valueOf(tem2[1]));                  
                }
            }
            System.out.println("Readed AABB" + "\n");
            while((temp3 = br3.readLine()) != null){
                String tem3[] = temp3.split("\t");
                if(tem3[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AABBDDset.add(Integer.valueOf(tem3[1]));                  
                }
            }
            System.out.println("Readed AABBDD" + "\n");
            //现在开始统计交并补集
            bw.write("Set111\tSet110\tSet101\tSet011\tSet001\tSet100\tSet010\n");
            int a1 = SetOperations_111(AAset,AABBset,AABBDDset);
            int a2 = SetOperations_110(AAset,AABBset,AABBDDset);
            int a3 = SetOperations_101(AAset,AABBset,AABBDDset);
            int a4 = SetOperations_011(AAset,AABBset,AABBDDset);
            int a5 = SetOperations_001(AAset,AABBset,AABBDDset);
            int a6 = SetOperations_100(AAset,AABBset,AABBDDset);
            int a7 = SetOperations_010(AAset,AABBset,AABBDDset);
            bw.write(a1 + "\t" + a2 + "\t" + a3 + "\t" + a4 + "\t" + a5 + "\t" + a6 + "\t" + a7 + "\n");
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    /** 这是对三个集合的操作,使用SetView完成对交集并集补集的操作
 * 集合的操作：交集、差集、并集  
 * Sets.intersection()交集
 * Sets.difference()差集
 * Sets.union()并集
 */
    
    public Integer SetOperations_111(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABB = Sets.intersection(AA, AABB);
        SetView<Integer> set111 = Sets.intersection(intersectionAA_AABB, AABBDD);
        return set111.size();
    }
    
    public Integer SetOperations_110(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABB = Sets.intersection(AA, AABB);
        SetView<Integer> set110 = Sets.difference(intersectionAA_AABB, AABBDD);
        return set110.size();
    }
    
    public Integer SetOperations_101(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABBDD = Sets.intersection(AA, AABBDD);
        SetView<Integer> set101 = Sets.difference(intersectionAA_AABBDD, AABB);
        return set101.size();
    }
    
    public Integer SetOperations_011(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAABB_AABBDD = Sets.intersection(AABB, AABBDD);
        SetView<Integer> set011 = Sets.difference(intersectionAABB_AABBDD, AA);
        return set011.size();
    }   
    
    public Integer SetOperations_001(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAA_AABB = Sets.union(AA, AABB);
        SetView<Integer> set001 = Sets.difference(AABBDD, unionAA_AABB);
        return set001.size();
    }  
    
    public Integer SetOperations_100(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAABB_AABBDD = Sets.union(AABB, AABBDD);
        SetView<Integer> set100 = Sets.difference(AA, unionAABB_AABBDD);
        return set100.size();
    }  
    
    public Integer SetOperations_010(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAA_AABBDD = Sets.union(AA, AABBDD);
        SetView<Integer> set010 = Sets.difference(AABB, unionAA_AABBDD);
        return set010.size();
    }  
    
      //在各个群体里面的分离情况进行统计
    //对Alineage来说，输入文件有3个，AA的个体的VCF文件，AABB的个体的VCF文件，AABBDD的个体的VCF文件
//    java -jar -Xms50g -Xmx70g /data2/xuebo/Projects/Speciation/javaCode/getVCFfile_info2.jar --file1 /data2/xuebo/Projects/Speciation/E2_all/macSS/chr${chr}.macSS.vcf.gz \
//	--file2 /data2/xuebo/Projects/Speciation/E2_all/macAABB/chr${chr}.macAABB.vcf.gz  --file3 /data2/xuebo/Projects/Speciation/E2_all/macAABBDD/chr${chr}.macAABBDD.vcf.gz \
//	--out /data2/xuebo/Projects/Speciation/E2_all/macSS/info2_chr${chr}.txt &
    //和上面的方法比较，不同的是，之前getVCFfile_info2得到的是数量，但是getVCFfile_info3得到的是位点
    
    public void getVCFfile_info3(String infileS1,String infileS2,String infileS3,String outfileS){
        try{
            //String temp = null;
            String temp1 = null;
            String temp2 = null;
            String temp3 = null;
            //int countall = 0;
            Set AAset = new HashSet();
            Set AABBset = new HashSet();
            Set AABBDDset = new HashSet();
            //HashMap<Integer, String> hashMap1 = new HashMap<Integer, String>();
            BufferedReader br1; //这里读的是AA个体的文件夹
            if (infileS1.endsWith("gz")) {
                br1 = IOUtils.getTextGzipReader(infileS1);
            } else {
                br1 = IOUtils.getTextReader(infileS1);
            }
            BufferedReader br2; //这里读的是AABB个体的文件夹
            if (infileS2.endsWith("gz")) {
                br2 = IOUtils.getTextGzipReader(infileS2);
            } else {
                br2 = IOUtils.getTextReader(infileS2);
            }
            BufferedReader br3; //这里读的是AABBDD个体的文件夹
            if (infileS3.endsWith("gz")) {
                br3 = IOUtils.getTextGzipReader(infileS3);
            } else {
                br3 = IOUtils.getTextReader(infileS3);
            }
            BufferedWriter bw = null; //写出文件，文件格式是txt格式         
            //开始读genome的VCF文件
            while((temp1 = br1.readLine()) != null){
                String tem1[] = temp1.split("\t");
                if(tem1[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AAset.add(Integer.valueOf(tem1[1]));                  
                }
            }
            System.out.println("Readed AA" + "\n");
            while((temp2 = br2.readLine()) != null){
                String tem2[] = temp2.split("\t");
                if(tem2[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AABBset.add(Integer.valueOf(tem2[1]));                  
                }
            }
            System.out.println("Readed AABB" + "\n");
            while((temp3 = br3.readLine()) != null){
                String tem3[] = temp3.split("\t");
                if(tem3[0].startsWith("#")){
                    //前边的不要,带有注释信息的那一行也不要
                }               
                else { //现在开始读的是没有注释的文件
                    AABBDDset.add(Integer.valueOf(tem3[1]));                  
                }
            }
            System.out.println("Readed AABBDD" + "\n");
            //现在开始统计交并补集
            //bw.write("Set111\tSet110\tSet101\tSet011\tSet001\tSet100\tSet010\n");
            SetView a1 = SetOperations2_111(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_111.txt");
                final List<Integer> list1 = new ArrayList<Integer>();  
                for(final Object value : a1){  
                    list1.add((Integer) value);  
                }  
                Collections.sort(list1);
                for(int n=0;n<list1.size();n++){
                    bw.write(list1.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            SetView a2 = SetOperations2_110(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_110.txt");
                final List<Integer> list2 = new ArrayList<Integer>();  
                for(final Object value : a2){  
                    list2.add((Integer) value);  
                }  
                Collections.sort(list2);
                for(int n=0;n<list2.size();n++){
                    bw.write(list2.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            SetView a3 = SetOperations2_101(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_101.txt");
                final List<Integer> list3 = new ArrayList<Integer>();  
                for(final Object value : a3){  
                    list3.add((Integer) value);  
                }  
                Collections.sort(list3);
                for(int n=0;n<list3.size();n++){
                    bw.write(list3.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            SetView a4 = SetOperations2_011(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_011.txt");
                final List<Integer> list4 = new ArrayList<Integer>();  
                for(final Object value : a4){  
                    list4.add((Integer) value);  
                }  
                Collections.sort(list4);
                for(int n=0;n<list4.size();n++){
                    bw.write(list4.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            SetView a5 = SetOperations2_001(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_001.txt");
                final List<Integer> list5 = new ArrayList<Integer>();  
                for(final Object value : a5){  
                    list5.add((Integer) value);  
                }  
                Collections.sort(list5);
                for(int n=0;n<list5.size();n++){
                    bw.write(list5.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            SetView a6 = SetOperations2_100(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_100.txt");
                final List<Integer> list6 = new ArrayList<Integer>();  
                for(final Object value : a6){  
                    list6.add((Integer) value);  
                }  
                Collections.sort(list6);
                for(int n=0;n<list6.size();n++){
                    bw.write(list6.get(n) + "\n");
                }  
            bw.flush();
            bw.close();
            Set a7 = SetOperations2_010(AAset,AABBset,AABBDDset);
            bw = IOUtils.getTextWriter(outfileS+"/SetOperations3_010.txt");
                final List<Integer> list7 = new ArrayList<Integer>();  
                for(final Object value : a7){  
                    list7.add((Integer) value);  
                }  
                Collections.sort(list7);
                for(int n=0;n<list7.size();n++){
                    bw.write(list7.get(n) + "\n");
                }  
            bw.flush();
            bw.close();

        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
     /** 这是对三个集合的操作,使用SetView完成对交集并集补集的操作,之前是为了得到交集，并集，补集的大小，现在是为了得到位点
 * 集合的操作：交集、差集、并集  
 * Sets.intersection()交集
 * Sets.difference()差集
 * Sets.union()并集
 */
    
    public SetView<Integer> SetOperations2_111(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABB = Sets.intersection(AA, AABB);
        SetView<Integer> set111 = Sets.intersection(intersectionAA_AABB, AABBDD);
        return set111;
    }
    
    public SetView<Integer> SetOperations2_110(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABB = Sets.intersection(AA, AABB);
        SetView<Integer> set110 = Sets.difference(intersectionAA_AABB, AABBDD);
        return set110;
    }
    
    public SetView<Integer> SetOperations2_101(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAA_AABBDD = Sets.intersection(AA, AABBDD);
        SetView<Integer> set101 = Sets.difference(intersectionAA_AABBDD, AABB);
        return set101;
    }
    
    public SetView<Integer> SetOperations2_011(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> intersectionAABB_AABBDD = Sets.intersection(AABB, AABBDD);
        SetView<Integer> set011 = Sets.difference(intersectionAABB_AABBDD, AA);
        return set011;
    }   
    
    public SetView<Integer> SetOperations2_001(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAA_AABB = Sets.union(AA, AABB);
        SetView<Integer> set001 = Sets.difference(AABBDD, unionAA_AABB);
        return set001;
    }  
    
    public SetView<Integer> SetOperations2_100(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAABB_AABBDD = Sets.union(AABB, AABBDD);
        SetView<Integer> set100 = Sets.difference(AA, unionAABB_AABBDD);
        return set100;
    }  
    
    public SetView<Integer> SetOperations2_010(Set<Integer> AA,Set<Integer> AABB,Set<Integer> AABBDD){
        SetView<Integer> unionAA_AABBDD = Sets.union(AA, AABBDD);
        SetView<Integer> set010 = Sets.difference(AABB, unionAA_AABBDD);
        return set010;
    }  
    
    //这个方法是给上面的代码产生的文件前面加上染色体
    public void getVCFfile_info4(String infileS,String outfileS){
         try{    
            String temp = null;
            BufferedReader siteFile = null;
            BufferedWriter bw = null;
            File f = new File(infileS);
            File[] fs = IOUtils.listRecursiveFiles(f);
            File[] sub = IOUtils.listFilesEndsWith(fs, ".txt");
            for(File fi:sub){
                siteFile = IOUtils.getTextReader(fi.toString());
                String Namelist[] = fi.toString().split("/");
                String chrName = Namelist[Namelist.length-2].split("_chr")[1];
                String fileName = Namelist[Namelist.length-1];
                System.out.print("It's " + chrName + "\n");
                bw = IOUtils.getTextWriter(outfileS + "/chr" + fileName);
                while((temp = siteFile.readLine()) != null){
                    bw.write(chrName + "\t" + temp + "\n");
                }   
                bw.flush();
                bw.close();
            } 
            
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    

   //这个方法的目的是去掉E2_all里面的空格
    
//    java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/cutALTblank.jar --file1 /data2/xuebo/Projects/Speciation/tree/barleybam/scan_barley/outchr${i}/VCF/chr${chr}.vcf \
//    --out /data2/xuebo/Projects/Speciation/tree/barleybam/scan_barley/outchr${i}/VCF/chr${chr}.vcf2  &
    
//     java -jar  -Xms20g -Xmx80g /data2/xuebo/Projects/Speciation/javaCode/cutALTblank2.jar --file1 /data2/xuebo/Projects/Speciation/E3/haveBlank/chr${chr}.all.vcf.gz \
//    --out /data2/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf  &
    public void cutALTblank(String infileS,String outfileS){
        try {
            String temp = null;
            BufferedReader br;
            if (infileS.endsWith("gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            } else {
                br = IOUtils.getTextReader(infileS);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //现在开始读VCF文件
            while((temp = br.readLine())!= null){
                if(temp.startsWith("#")){
                    bw.write(temp + "\n");
                }else{
                   String tem[] = temp.split("\t");
                   String barleyhaplo = tem[tem.length-1].split(":")[0];
                   if(barleyhaplo.equals("./.")){
                       
                   }else{
                       bw.write(temp + "\n");
                   }                             
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
