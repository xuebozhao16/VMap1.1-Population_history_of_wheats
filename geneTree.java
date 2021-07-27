/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import static EvolutionWheat.ForVcftoolsGroup.getTextReader;
import static EvolutionWheat.ForVcftoolsGroup.listFilesEndsWith;
import static EvolutionWheat.ForVcftoolsGroup.listRecursiveFiles;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class geneTree {
    public geneTree(String infileS1,String infileS2,String outfileS){
        //this.C20_countgeneSNP(infileS1, infileS2, outfileS);
        //this.C21_getgeneSNP(infileS1, infileS2, outfileS);
        //this.C22_fasta2phy(infileS1,infileS2,outfileS);
        //this.C24_addlinenum(infileS1, infileS2, outfileS);
        //this.C25_subspeciesname(infileS1, infileS2, outfileS);
        //this.C28_musclephy2raxmlphy(infileS1, infileS2, outfileS);
        this.C31_getTreetopologyNUM(infileS1, infileS2, outfileS);
    }
    public geneTree(String infileS,String outfileS){
        //this.C23_fasta2phy(infileS, outfileS);
        //this.C29_forgene10Mtree_head(infileS, outfileS);
        //this.C30_getTreetopology(infileS, outfileS);
        this.C32_RAxML_bestTree(infileS, outfileS);
    }
    
    public geneTree(String infileS1,String infileS2,String outfileS1,String outfileS2,String outfileS3){
        this.get111_bedfile(infileS1, infileS2, outfileS1, outfileS2, outfileS3);
    }
    
    //现在这个方法是统计每个基因里面有多少个SNP //染色体已经分好,这里不用再去判断染色体 //重叠的gene算，但是包含的基因不算
//        String infileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/V1_1/GeneLulab1_1onlyGene/chr36_GeneLulab1_1.txt";
//        String infileS2 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/test/chr36.allgenic.vcf";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/geneTree/test/chr36_Genecount.txt";
//        new geneTree(infileS1,infileS2,outfileS);
    // 输出文件格式 36	22460	24375	TraesCS6D02G355900	8
    public void C20_countgeneSNP(String infileS1,String infileS2,String outfileS){
        try {
            String temp = null;
            String temp2 = null;
            BufferedReader br1;
            BufferedReader br2;
            int i = 0;
            int count = 0;
            String chr = null;
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
            List<Integer> Pos1List = new ArrayList <>();
            List<Integer> Pos2List = new ArrayList <>();
            List<String> geneIDList = new ArrayList <>();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //现在开始读gene注释文件，//文件格式是 1	40098	70338	TraesCS1A02G000100
            int pos1 = 0;
            int pos2 = 0;
            while((temp = br1.readLine())!= null){
                String tem[] = temp.split("\t");
                //现在要判断重叠和包含的基因
                if(Integer.valueOf(tem[1])<= pos2 && Integer.valueOf(tem[2])>=pos2){//这个是交叉的情况
                    Pos1List.add(Integer.valueOf(tem[1]));
                    Pos2List.add(Integer.valueOf(tem[2]));
                    geneIDList.add(tem[3]);
                    pos1 = Integer.valueOf(tem[1]); //输出之后要赋值
                    pos2 = Integer.valueOf(tem[2]);
                }
                else if(Integer.valueOf(tem[1])<= pos2 && Integer.valueOf(tem[2])<=pos2) continue;//这个是上一个包含这个
                else if(Integer.valueOf(tem[1])> pos2) {
                    Pos1List.add(Integer.valueOf(tem[1]));
                    Pos2List.add(Integer.valueOf(tem[2]));
                    geneIDList.add(tem[3]);
                    pos1 = Integer.valueOf(tem[1]); //输出之后要赋值
                    pos2 = Integer.valueOf(tem[2]);
                }
            }
            System.out.println(Pos1List.size() + "\n");
            //System.out.println(Pos2List.size() + "\n");
            //System.out.println(geneIDList.size() + "\n");
            //现在开始读VCF文件
            while((temp2 = br2.readLine())!= null){
                //String tem2[] = temp2.split("\t");
                if(!temp2.startsWith("#")){
                    chr = temp2.split("\t")[0];
                    int pos = Integer.valueOf(temp2.split("\t")[1]);
                    //现在开始做判断
                    //if(pos < Pos1List.get(i)) continue;
                    if (pos >= Pos1List.get(i) && pos <= Pos2List.get(i)){
                        count = count + 1;
                    }
                    else if (pos >= Pos1List.get(i+1) && pos <= Pos2List.get(i+1)){
                        bw.write(chr + "\t" + Pos1List.get(i) + "\t" + Pos2List.get(i) + "\t" + geneIDList.get(i) + "\t" + count + "\n");
                        count = 1; //因为这个时候已经是下一个基因了
                        i = i + 1;
                    }
                    else if (pos > Pos2List.get(i+1)){
                        bw.write(chr + "\t" + Pos1List.get(i) + "\t" + Pos2List.get(i) + "\t" + geneIDList.get(i) + "\t" + "0" + "\n");
                        i = i + 1;
                    }
                }       
            }
            bw.write(chr + "\t" + Pos1List.get(i) + "\t" + Pos2List.get(i) + "\t" + geneIDList.get(i) + "\t" + count + "\n");
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //这个方法是为了得到几万个gene的SNP文件
    
    // 输入文件格式 36	22460	24375	TraesCS6D02G355900	8
    // 输入文件格式 36	22460	24375	3452
    public void C21_getgeneSNP(String infileS1,String infileS2,String outfileS){
        try {
            String temp = null;
            String temp2 = null;
            BufferedReader br1;
            BufferedReader br2;
            int i = 0;
            int nnn = 0;
            //StringBuilder pergeneSNP = null;
            String headline = null;
            String geneID = null;
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
            List<Integer> Pos1List = new ArrayList <>();
            List<Integer> Pos2List = new ArrayList <>();
            List<String> geneIDList = new ArrayList <>();
            List<String> geneSNP = new ArrayList <>();
            BufferedWriter bw = null;
            //现在开始读gene注释文件，//文件格式是 1	40098	70338	TraesCS1A02G000100  8
            while((temp = br1.readLine())!= null){
                nnn = nnn+1;
                String tem[] = temp.split("\t");
                Pos1List.add(Integer.valueOf(tem[1]));
                Pos2List.add(Integer.valueOf(tem[2]));
                geneIDList.add("C"+ tem[0] + "gene" + nnn);
            }
            System.out.println(Pos1List.size() + "\n");
            //System.out.println(Pos2List.size() + "\n");
            //System.out.println(geneIDList.size() + "\n");
            //现在开始读VCF文件
            while((temp2 = br2.readLine())!= null){
                //String tem2[] = temp2.split("\t");
                if(temp2.startsWith("#C")){
                    headline = temp2;
                }
                if(!temp2.startsWith("#")){
                    if(i<Pos2List.size()-1){
                        //chr = temp2.split("\t")[0];
                        int pos = Integer.valueOf(temp2.split("\t")[1]);
                        //现在开始做判断
                        //if(pos < Pos1List.get(i)) continue;
                        if (pos >= Pos1List.get(i) && pos <= Pos2List.get(i)){
                            //System.out.println(pergeneSNP);
                            geneSNP.add(temp2);
                        }
                        else if(pos > Pos2List.get(i) && pos < Pos1List.get(i+1)) continue;
                        else if (pos >= Pos1List.get(i+1) && pos <= Pos2List.get(i+1)) {
                            //因为这个时候SNP已经是下一个Pos了
                            geneID = geneIDList.get(i) ;
                            //bw = IOUtils.getTextWriter(outfileS+"/"+ geneID +".vcf");
                            if(geneSNP.size() > 50){
                                bw = IOUtils.getTextGzipWriter(outfileS+"/"+ geneID +".vcf.gz"); //这是输出的就是压缩文件
                                bw.write(headline + "\n");
                                for(int j=0;j<geneSNP.size();j++){
                                    bw.write(geneSNP.get(j) + "\n");
                                }
                                bw.flush();
                                bw.close();
                            }
                            //现在要清空geneSNP啦  List<String> geneSNP = new ArrayList <>();
                            geneSNP.clear();
                            geneSNP.add(temp2);
                            i = i + 1;
                        }
                        else if (pos > Pos2List.get(i+1)){
                            i = i + 1;
                        }
                    }
                    
                }   
            }
            geneID = geneIDList.get(i) ;
            if(geneSNP.size() > 50){
            bw = IOUtils.getTextGzipWriter(outfileS+"/"+ geneID +".vcf.gz");
                bw.write(headline + "\n");
                for(int j=0;j<geneSNP.size();j++){
                    bw.write(geneSNP.get(j) + "\n");
                }
                bw.flush();
                bw.close();
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //这个方法是从fasta文件到phy文件,这个方法是已知个体的数量
    public void C22_fasta2phy(String infileS,String samplenum,String outfileS){
        try {
            String temp = null;
            String sample = null;
            BufferedReader br = IOUtils.getTextReader(infileS); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            int i = 3;
            //前两行单独弄
            String line1 = br.readLine().replace(">", "");
            String line2 = br.readLine();
            bw.write(samplenum + "\t" + line2.length() + "\n");
            bw.write(line1 + "\t" + line2 + "\n");
            while((temp = br.readLine())!= null){
               if(i%2 == 1){
                   sample = temp.replace(">", "");
               }else{
                   bw.write(sample + "\t" + temp + "\n");
               }
               i = i+1;
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //这个方法是从fasta文件到phy文件,这个方法是不知道个体的数量
    public void C23_fasta2phy(String infileS,String outfileS){
        try {
            String temp = null;
            //String line = null;
            String sample = null;
            int countsample = 0;
            int templength = 0;
            BufferedReader br = IOUtils.getTextReader(infileS); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);       
            int i = 1;
            List<String> fastaList = new ArrayList <>();       
            while((temp = br.readLine())!= null){
               if(i%2 == 1){
                    sample = temp.replace(">", "");
                    countsample = countsample + 1;
               }else{
                    templength = temp.length();
                    fastaList.add(sample + "\t" + temp + "\n");
               }
               i = i + 1;
            }
            bw.write(countsample + "\t" + templength + "\n");
            for(int j = 0;j<fastaList.size();j++){
                bw.write(fastaList.get(j));
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    
    //这个方法是加上行数，
     public void C24_addlinenum(String infileS1,String infileS2,String outfileS){
        try {
            String temp = null;
            String temp2 = null;
            int linenum = 0;
            BufferedReader br1 = IOUtils.getTextReader(infileS1); 
            BufferedReader br2 = IOUtils.getTextReader(infileS2); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //HashMap<String, String> hashMap1 = new HashMap<String, String>();
            List<String> num = new ArrayList <>();
            List<String> name = new ArrayList <>();
            while((temp = br1.readLine())!= null){
                String tem[] = temp.split("\t");
                //hashMap1.put(tem[1],tem[0]);
                num.add(tem[0]);
                name.add(tem[1]);
            }
            while((temp2 = br2.readLine())!= null){
                String newtemp2 = temp2;
                for(int i = 0;i< name.size();i++){
                    String value = name.get(i);
                    String newvalue = num.get(i);
                    //System.out.println(value + "\n");
                    newtemp2 = newtemp2.replace(value, newvalue);
                    //newtemp2 = temp2.replace(value, Integer.toString(i));
                }
                bw.write("tree STATE_" + linenum*1000 + " = " + newtemp2 + "\n");
                linenum = linenum + 1;
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
     
     
     //这个方法是得到不同的subspecies的个体，出来的是27个subspecies的个体的名字
     //输入文件是 1:Wild_einkorn Domesticated_einkorn ~   2: A001	Wild_einkorn    A002	Wild_einkorn
     //输出文件是 Wild_einkorn.txt ~ ~
     public void C25_subspeciesname(String infileS1,String infileS2,String outfileS){
        try {
            String temp = null;
            String temp2 = null;
            int count = 0;
            BufferedReader br1 = IOUtils.getTextReader(infileS1); 
            BufferedReader br2 = IOUtils.getTextReader(infileS2); 
            BufferedWriter bw = null;
            List<String> subList = new ArrayList <>(); 
            List<String> sampleList = new ArrayList <>(); 
            //现在开始一行行的读所有个体的文件，放到Array里面
            while((temp = br1.readLine())!= null){
                String tem[] = temp.split("\t");
                subList.add(tem[1]);
                sampleList.add(tem[0]);
            }
            // 读的是有哪些subspecies的文件
            while((temp2 = br2.readLine())!= null){
                int subcount = 0;
                bw = IOUtils.getTextWriter(outfileS+"/sub_"+ temp2 +".txt");
                for(int j = 0;j<subList.size();j++){
                    if(temp2.equals(subList.get(j))){
                        count = count + 1;
                        subcount = subcount + 1;
                        bw.write(sampleList.get(j) + "\n");
                    }
                }
                //System.out.println(temp2 + "\" + subcount + ");
                bw.flush();
                bw.close();
            }
            System.out.println(count);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
     
     //之前的方法行不通，现在开始从头组装/data2/xuebo/Projects/Speciation/geneTree/assembly
     //得到 1:1:1的基因的bed文件,上游150bp ,下游150bp
     public void get111_bedfile(String infileS1,String infileS2,String outfileS1,String outfileS2,String outfileS3){
        try {
            String temp = null;
            String temp2 = null;
            int count = 0;
            BufferedReader br1 = IOUtils.getTextReader(infileS1); 
            BufferedReader br2 = IOUtils.getTextReader(infileS2); 
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1); 
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2); 
            BufferedWriter bw3 = IOUtils.getTextWriter(outfileS3); 
            Set A = new HashSet();
            Set B = new HashSet();
            Set D = new HashSet();
            //1：1：1的文件
            while((temp = br1.readLine())!= null){
                String tem[] = temp.split("\t");
                A.add(tem[0]);
                B.add(tem[1]);
                D.add(tem[2]);
            }
            // 读的是 42	TraesCS7D02G515400	162055927	162058743	-  文件
            while((temp2 = br2.readLine())!= null){
                String tem[] = temp2.split("\t");
                int pos1 = Integer.valueOf(tem[2])-150;
                int pos2 = Integer.valueOf(tem[3])+150;
                if(!A.add(tem[1])){
                    //bw1.write(tem[0] + "\t" + tem[2] + "\t" + tem[3] + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                    bw1.write(tem[0] + "\t" + pos1 + "\t" + pos2 + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                }else{
                    A.remove(tem[1]);
                }
                if(!B.add(tem[1])){
                    //bw2.write(tem[0] + "\t" + tem[2] + "\t" + tem[3] + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                    bw2.write(tem[0] + "\t" + pos1 + "\t" + pos2 + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                }else{
                    B.remove(tem[1]);
                }
                if(!D.add(tem[1])){
                    //bw3.write(tem[0] + "\t" + tem[2] + "\t" + tem[3] + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                    bw3.write(tem[0] + "\t" + pos1 + "\t" + pos2 + "\t" + tem[1] + "\t" + tem[4] + "\n" );
                }else{
                    D.remove(tem[1]);
                } 
            }
            bw1.flush();
            bw2.flush();
            bw3.flush();
            bw1.close();
            bw2.close();
            bw3.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
     //这个方法是从比对软件中得到的phy文件到raxml可以识别的phy文件,
     //这个方法是知道个体的数量,但是名字不是很全，要自己补齐
    public void C28_musclephy2raxmlphy(String infileS1,String infileS2,String outfileS){
        try {
            String temp1 = null;
            String temp = null;
            StringBuilder line  = new StringBuilder();
            BufferedReader br1 = IOUtils.getTextReader(infileS1); 
            BufferedReader br2 = IOUtils.getTextReader(infileS2);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            HashMap<String,String> hashMap1 = new HashMap<String,String>();
            Set ABD = new HashSet();
            List<String> phyList = new ArrayList <>();
            while((temp1 = br1.readLine())!= null){
                hashMap1.put(temp1.split("\t")[0], temp1.split("\t")[1]);
                ABD.add(temp1.split("\t")[0]);
            }
            bw.write(br2.readLine() + "\n");
            while((temp = br2.readLine())!= null){
                String tem[] = temp.split(" ");
                String head = null;
                if(tem[0].contains(" ")){
                    head = tem[0].split(" ")[0];   
                }else{
                    head = tem[0];
                }           
                for(int i = 0;i<tem.length;i++){
                    if(i==0){
                        if(!ABD.add(head)){
                            phyList.add(line.toString());
                            line.delete( 0, line.length() ); 
                            phyList.add(hashMap1.get(head));
                            //System.out.println(line+"\n");
                            
                            
                        }else{
                            line.append(tem[0]);
                            ABD.remove(tem[0]);
                        }
                    }else{
                        line.append(tem[i]);
                    } 
                }
            }
            phyList.add(line.toString());
            int j = 0;
            for(j = 0;j<phyList.size();j++){
                if(phyList.get(j).length()!=0){
                    if(j%2==0){
                        bw.write(phyList.get(j) + "\n");
                    }else{
                        bw.write(phyList.get(j) + "\t");
                    }
                }   
            }       
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //读文件夹里面的所有树，之后吧所有的树前面加上这是哪个染色体哪个10M的数据
    public void C29_forgene10Mtree_head(String infileS,String outfileS){
         try{    
            String temp = null; 
            String temporder = null;
            BufferedReader xpclrFile = null;
            //BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //HashMap<String, Double> hashMapMean = new HashMap<String, Double>();
            //bw.write("taxa" + "\t" + "pi" + "\n");
            File f = new File(infileS);
            File[] fs = listRecursiveFiles(f);
            File[] sub = listFilesEndsWith(fs, ".tree");
            for(File fi:sub){
                int i = 0;
                double sum = 0;
                xpclrFile = getTextReader(fi.toString());
                String taxaNamelist[] = fi.toString().split("/");
                //String taxaName1 = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
                String taxaName1 = taxaNamelist[taxaNamelist.length-1].split("\\.")[1];
                String chr = taxaName1.split("_")[0];
                int pos1 = Integer.valueOf(taxaName1.split("_")[1]);
                int pos2 = pos1 + 10000000;
                System.out.print("It's " + taxaName1 + "\n");
                while((temp = xpclrFile.readLine()) != null){
                    bw.write(chr + "\t" + pos1 + "\t" + pos2 + "\t" + temp + "\n" );
                }
            } 
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    //geneTree_10M,把里面的拓扑结构提取出来
    public void C30_getTreetopology(String infileS,String outfileS){
        try {
            String temp = null;
            //int count = 0;
            BufferedReader br = IOUtils.getTextReader(infileS); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder line  = new StringBuilder();
            while((temp = br.readLine())!= null){
                String tem[] = temp.split("\t");
                String topology[] = tem[3].split(",");
                for(int i=0;i<topology.length;i++){
                    line.append(topology[i].split(":")[0] + ":");
                }
                String  line1 = line.toString();
                String  line2 = line1.replace("(", "");
                String line3 = line2.substring(0,line2.length()-1);
                bw.write(tem[0] + "_" + tem[1] + "_" + tem[2] + "\t" + line3 + "\n");
                line.delete(0, line.length());
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //把上面的亚种代替为数字
    public void C31_getTreetopologyNUM(String infileS1,String infileS2,String outfileS){
        try {
            String temp = null;
            String temp1 = null;
            BufferedReader br1 = IOUtils.getTextReader(infileS1); 
            BufferedReader br2 = IOUtils.getTextReader(infileS2); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder line  = new StringBuilder();
            HashMap<String, String> hashMap1 = new HashMap<String, String>();
            Set A = new HashSet();
            while((temp1 = br1.readLine())!= null){
                hashMap1.put(temp1.split("\t")[1], temp1.split("\t")[0]);
                A.add(temp1.split("\t")[1]);
            }
            while((temp = br2.readLine())!= null){
                String tem[] = temp.split("\t");
                String topology[] = tem[1].split(":");
                for(int i=0;i<topology.length;i++){
                    if(!A.add(topology[i])){
                        line.append(hashMap1.get(topology[i]) + "\t");
                    }else{
                        A.remove(topology[i]);
                    }
                }
                String  line1 = line.toString();    
                String line2 = line1.substring(0,line1.length()-1);
                bw.write(tem[0] + "\t" + line2 + "\n");
                line.delete(0, line.length());
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //这个方法是对文件夹里面RAxML_bestTree开头的文件进行操作，
    //例如 文件夹RAxML_bestTree.35_90000000.tree，输出的文件是35\t8000000\t90000000\t文件夹里面的内容
    public void C32_RAxML_bestTree(String infileS,String outfileS){
         try{    
            String temp = null; 
            BufferedReader piFile = null;
            //BufferedWriter bw = null;
            File f = new File(infileS);
            File[] fs = listRecursiveFiles(f);
            File[] sub = listFilesEndsWith(fs, ".tree");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //Set<String> features  = new HashSet();
            //Map<String,ArrayList<String[]>> data = new HashMap();
            for(File fi:sub){
                piFile = getTextReader(fi.toString());
                String taxaNamelist[] = fi.toString().split("/");
                String taxaName = taxaNamelist[taxaNamelist.length-1].split("\\.")[1];
                System.out.print("It's " + taxaName + "\n");
                while((temp = piFile.readLine()) != null){
                    //System.out.print(taxaName.split("_")[1] + "\n");
                    String pos1 = taxaName.split("_")[1];
                    String chr  = taxaName.split("_")[0];
                    int pos2 = Integer.valueOf(taxaName.split("_")[1]) + 10000000;
                    String tree  = temp;
                    bw.write(chr+ "\t" + pos1 + "\t" + pos2 + "\t" + tree + "\n");
                }
            } 
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    
    
    
}
