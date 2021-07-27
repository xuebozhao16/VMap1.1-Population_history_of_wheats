/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */ //这个代码是为了处理画树的代码
public class tree {
    public tree(String outfile){
        //this.gethapScan_barley_taxaBAMfile(outfile);
        //this.gethapScan_Ae_taxaBAMfile(outfile);
        this.gethapScan_Tu_taxaBAMfile(outfile);
    }
    
    public tree(String infileS,String outfileS){
        //this.getwithBarleyVCF(infileS, outfileS);
        //this.cutALTblank(infileS, outfileS);
        this.C26_getwithBarley_Ae_VCF(infileS, outfileS);
    }
    
    public tree(String infileS1,String infileS2,String outfileS){
        this.C17_getAABB_AABBDD_segregateVCF(infileS1, infileS2, outfileS);
    }
    public tree(String infileS1,String infileS2,String infileS3,String outfileS){
        this.C18_getAABB_AABBDD_segregateVCF(infileS1, infileS2, infileS3, outfileS);
    }
    
    
    //这个代码是为了获得barley做outgroup，把barley做scan,得到barley的VCF文件，之后使用BCFtools合并起来
    //barley的bam文件是分染色体的
//        String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatScan/hapScanner2/Parameters_hapScanner2_chr1.txt";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/tree/parameters_file";
//        new syntenic_Sgenome(infileS,outfileS);
    
    public void gethapScan_barley_taxaBAMfile(String outfileS){
        try {            
            BufferedWriter bw = null;
            String[] te = null;          
            for(int i = 1; i < 43;i++ ){                
                bw = IOUtils.getTextWriter(outfileS+"/TaxaRefBamChr"+i);
                bw.write("Taxa" + "\t" + "Reference" + "\t" + "BamPath" + "\n");
                bw.write("barley" + "\t" + "/data1/home/xuebo/Project/reference/abd_iwgscV1.fa.gz" + "\t" + 
                                "/data2/xuebo/Projects/Speciation/tree/barleybam/chr"+i+".sorted.bam" + "\n");
                bw.flush();
                bw.close();
            }
            
        } catch (Exception e){
            e.printStackTrace();
        }
    }
    
    public void gethapScan_Ae_taxaBAMfile(String outfileS){
        try {            
            BufferedWriter bw = null;
            String[] te = null;    
            int chrA[] = {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40};
            for(int i = 0; i < chrA.length;i++ ){                
                bw = IOUtils.getTextWriter(outfileS+"/TaxaRefBamChr"+chrA[i]);
                bw.write("Taxa" + "\t" + "Reference" + "\t" + "BamPath" + "\n");
                bw.write("Ae" + "\t" + "/data1/home/xuebo/Project/reference/abd_iwgscV1.fa.gz" + "\t" + 
                                "/data2/xuebo/Projects/Speciation/introgression/outgroup/Ae/chr"+chrA[i]+".sorted.bam" + "\n");
                bw.flush();
                bw.close();
            }
            
        } catch (Exception e){
            e.printStackTrace();
        }
    }
    
    public void gethapScan_Tu_taxaBAMfile(String outfileS){
        try {            
            BufferedWriter bw = null;
            String[] te = null;    
            int chrA[] = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};
            for(int i = 0; i < chrA.length;i++ ){                
                bw = IOUtils.getTextWriter(outfileS+"/TaxaRefBamChr"+chrA[i]);
                bw.write("Taxa" + "\t" + "Reference" + "\t" + "BamPath" + "\n");
                bw.write("Tu" + "\t" + "/data1/home/xuebo/Project/reference/abd_iwgscV1.fa.gz" + "\t" + 
                                "/data2/xuebo/Projects/Speciation/introgression/outgroup/Tu/chr"+chrA[i]+".sorted.bam" + "\n");
                bw.flush();
                bw.close();
            }
            
        } catch (Exception e){
            e.printStackTrace();
        }
    }
    
    //这个方法是为了row_chr.withBarley.vcf里面barley是./.的去掉
    public void getwithBarleyVCF(String infileS,String outfileS){
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
    
    
    //这个方法是为了row_chr.withBarleyAe.vcf里面barley和Ae是./.的去掉
    public void C26_getwithBarley_Ae_VCF(String infileS,String outfileS){
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
                   String Aehaplo = tem[tem.length-2].split(":")[0];
                   if(!barleyhaplo.equals("./.") &&  !Aehaplo.equals("./.") ){
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
    
    //这个方法的目的是去掉/data2/xuebo/Projects/Speciation/tree/barleybam/scan_barley_vcf 里面的一个空格，这里面原来是两个空格的
    ///data2/xuebo/Projects/Speciation/E2_all/chr${chr}.all.vcf.gz 因为这个里面是有一个空格的，这样都有一个空格才可以merge起来
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
                   //String cutALTblank = tem[4].substring(0, 2); //这个方法的目的是去掉/data2/xuebo/Projects/Speciation/tree/barleybam/scan_barley_vcf 里面的一个空格，这里面原来是两个空格的  cutALTblank.jar
                   String cutALTblank = tem[4].substring(0, 1); //这个方法的目的是去掉/data2/xuebo/Projects/Speciation/E3/haveBlank 里面的空格，得到标准的E3文件 cutALTblank2.jar
                   temp = temp.replace(tem[4], cutALTblank);
                   bw.write(temp + "\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //这个方法是为了得到AABB和AABBDD中的所有分离位点和总的位点的一部分得到一个Set，这里面全是AABB和AABBDD的分离的位点
    public void C17_getAABB_AABBDD_segregateVCF(String infileS1,String infileS2,String outfileS){
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
            Set Set1 = new HashSet();
            //List<String> list2 = new ArrayList<String> ();  
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //AABB和AABBDD中的所有分离位点的VCF文件,把里面的所有位点都放进Set里面
            while((temp1 = br1.readLine())!= null){
                if(!temp1.startsWith("#")){
                    Set1.add(temp1.split("\t")[1]);
                }
            } 
            System.out.println(Set1.size() + "\n");
            //现在Set1里面就是AABB和AABBDD中的所有分离位点
            //现在挑选Set1里面的SNP
            while((temp2 = br2.readLine())!= null){
                if(temp2.startsWith("#")){
                    bw.write(temp2 + "\n");
                }else{
                    String tem[] = temp2.split("\t");
                    if(!Set1.add(tem[1])){
                        bw.write(temp2 + "\n");
                    }else{
                        Set1.remove(tem[1]);
                    }           
                }
            }   
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //这个方法是为了得到AABB和AABBDD中的所有分离位点和总的位点的一部分得到一个Set，这里面AABB和AABBDD的分离的位点比较多
    public void C18_getAABB_AABBDD_segregateVCF(String infileS1,String infileS2,String infileS3,String outfileS){
        try {
            String temp1 = null;
            String temp2 = null;
            String temp3 = null;
            BufferedReader br1;
            BufferedReader br2;
            BufferedReader br3;
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
            if (infileS3.endsWith("gz")) {
                br3 = IOUtils.getTextGzipReader(infileS3);
            } else {
                br3 = IOUtils.getTextReader(infileS3);
            }
            Set Set1 = new HashSet();
            Set Set2 = new HashSet();
            Set Set3 = new HashSet();
            List<String> list2 = new ArrayList<String> ();  
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //AABB和AABBDD中的所有分离位点的VCF文件,把里面的所有位点都放进Set里面
            while((temp1 = br1.readLine())!= null){
                if(!temp1.startsWith("#")){
                    Set1.add(temp1.split("\t")[1]);
                }
            } 
            System.out.println(Set1.size() + "\n");
            while((temp2 = br2.readLine())!= null){
                if(!temp2.startsWith("#")){
                    Set2.add(temp2.split("\t")[1]);
                }
            }
            System.out.println(Set2.size() + "\n");
            list2.addAll(Set2);  
            System.out.println(list2.get(0) + "\n");
            for(int i = 0; i < list2.size()/6; i++){ 
                //String AA = list2.get(i);
                Set1.add(list2.get(i));
            }
            System.out.println(Set1.size() + "\n");
            //现在Set1里面就是AABB和AABBDD中的所有分离位点和总的位点的一部分
            //现在挑选Set1里面的SNP
            while((temp3 = br3.readLine())!= null){
                if(temp3.startsWith("#")){
                    bw.write(temp3 + "\n");
                }else{
                    String tem[] = temp3.split("\t");
                    if(!Set1.add(tem[1])){
                        bw.write(temp3 + "\n");
                    }else{
                        Set1.remove(tem[1]);
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
