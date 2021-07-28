package speciation;

import static EvolutionWheat.ForVcftoolsGroup.*;
import static com.sun.scenario.Settings.set;
//import format.table.RowTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
//import utils.IOUtils;

import java.io.*;
import java.lang.reflect.Array;
import java.util.HashSet;
import java.util.Set;

import static java.lang.reflect.Array.*;
import java.util.ArrayList;
import static java.util.Arrays.sort;
import java.util.HashMap;
import java.util.List;

public class split_time {
    public split_time(String infileS,String outfileS){
          //this.get_singlecopy_gene(infileS, outfileS);
          this.get_singlecopy_gene2(infileS, outfileS);
    }
    public split_time(String infileS1,String infileS2,String outfileS1,String outfileS2,String outfileS3){
        this.get_singlecopy_gene_bedfile(infileS1, infileS2, outfileS1, outfileS2, outfileS3);
    }
    public split_time(String infileS1,String infileS2,String outfileS){
        //this.C39_get_xmlfile_forbeast(infileS1, infileS2, outfileS);
        //this.FortheRealPos_gff(infileS1, infileS2, outfileS);
        //this.get_singlecopy_gene3_more1000(infileS1, infileS2, outfileS);
        //this.C42_move1_1genefile(infileS1, infileS2, outfileS);
        //this.C39_get_xmlfile_forbeast2(infileS1, infileS2, outfileS);
        //this.C39_get_xmlfile_forbeast3_A(infileS1, infileS2, outfileS);
        this.C39_get_xmlfile_forbeast3_B(infileS1, infileS2, outfileS);
        //this.C39_get_xmlfile_forbeast3_D(infileS1, infileS2, outfileS);
    }

    //这个方法是对之前Swiftorth的结果进行整理，找出single copy的基因
//    String infileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene.txt";
//    String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.txt";
//    new split_time(infileS,outfileS);
    public void get_singlecopy_gene(String infileS,String outfileS){
        try {
            BufferedWriter bw = null;
            String temp = "";
            String[] te = null;
            BufferedReader br = IOUtils.getTextReader(infileS);
            bw = IOUtils.getTextWriter(outfileS);
            //bw.write("Hv\tAe\tTu\tTdA\tTdb\tTtA\tTtB\tTaA\tTaB\tTaD");
            //bw.newLine();
            while((temp = br.readLine())!=null){
                Set<String> geneset = new HashSet<String>();
                te = temp.split("\t");
                for(int i = 0;i<te.length;i++){
                    if(te[i].contains(".")){
                        geneset.add(te[i].split("\\.")[0]);
                    }else{
                        geneset.add(te[i]);
                    }
                }
                if (geneset.size() == 10){
                    StringBuilder temmp = new StringBuilder();
                    String[] geneArray = new String[geneset.size()];
                    geneset.toArray(geneArray);
                   sort(geneArray);
                   for(int j = 0;j<geneArray.length;j++){
                       temmp.append(geneArray[j] + "\t");
                   }
                    bw.write(temmp.substring(0,temmp.length()-1) + "\n");
                }
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    //这个方法是得到single copy的基因的bed file
    //这里面的输入文件是两个，一个是含有十列的single copy gene list，一个是小麦的gff3文件
    //输出的是一个文件夹，里面包含有三个文件，分别是A，B，D的bed文件。Chr pos1 pos2 geneID
//    String infileS1 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100.txt";
//        String infileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/V1_1/GeneLulab1_1onlyGene.txt";
//        String outfileS1 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100A.bed";
//        String outfileS2 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100B.bed";
//        String outfileS3 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100D.bed";
//        new split_time(infileS1,infileS2,outfileS1,outfileS2,outfileS3);
    public void get_singlecopy_gene_bedfile(String infileS1,String infileS2,String outfileS1,String outfileS2,String outfileS3){
        try {
            String temp = null;
            String temp2 = null;
            String[] te = null;
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedReader br2 = IOUtils.getTextReader(infileS2);
            //bw.write("Hv\tAe\tTu\tTdA\tTdb\tTtA\tTtB\tTaA\tTaB\tTaD");
            //bw.newLine();
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            BufferedWriter bw3 = IOUtils.getTextWriter(outfileS3);
            Set<String> geneA = new HashSet<String>();
            Set<String> geneB = new HashSet<String>();
            Set<String> geneD = new HashSet<String>();
            while((temp = br.readLine())!=null){
                te = temp.split("\t");
                for(int i = 0;i<te.length;i++) {
                    if (te[i].contains("TaA")) {
                        //System.out.println(te[i]);
                        geneA.add(te[i].split("\\|")[1]);
                    }
                    if (te[i].contains("TaB")) {
                        geneB.add(te[i].split("\\|")[1]);
                    }
                    if (te[i].contains("TaD")) {
                        geneD.add(te[i].split("\\|")[1]);
                    }
                }
            }
            while((temp2 = br2.readLine())!=null){
                 String[] te2 = temp2.split("\t");
                 if(!geneA.add(te2[1])){
                     bw1.write(te2[0] + "\t" + (Integer.valueOf(te2[2])-1)+ "\t" +  te2[3] + "\t" +  te2[4]  + "\t" +  te2[1] + "\n");
                 }else{
                     geneA.remove(te2[1]);
                 }
                 if (!geneB.add(te2[1])){
                     bw2.write(te2[0] + "\t" + (Integer.valueOf(te2[2])-1)+ "\t" +  te2[3] + "\t" +  te2[4]  + "\t" +  te2[1] + "\n");
                 }else{
                     geneB.remove(te2[1]);
                 }
                 if (!geneD.add(te2[1])){
                     bw3.write(te2[0] + "\t" + (Integer.valueOf(te2[2])-1)+ "\t" +  te2[3] + "\t" +  te2[4]  + "\t" +  te2[1] + "\n");
                 }else{
                     geneD.remove(te2[1]);
                }
            }
            bw1.flush();
            bw2.flush();
            bw3.flush();
            bw1.close();
            bw2.close();
            bw3.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    //现在的目的是为了得到100个基因的bed文件，一行一个
//            String infileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100D.bed";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.gene100D";
//        new split_time(infileS,outfileS);
    public void get_singlecopy_gene2(String infileS,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp = "";
            String[] te = null;
            int line = 1;
            BufferedReader br = IOUtils.getTextReader(infileS);
            //bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                //BufferedWriter bw = IOUtils.getTextWriter(outfileS+"/singlecopygene_A_"+ line+".bed");
                //BufferedWriter bw = IOUtils.getTextWriter(outfileS+"/singlecopygene_B_"+ line+".bed");
                BufferedWriter bw = IOUtils.getTextWriter(outfileS+"/singlecopygene_D_"+ line+".bed");
                bw.write(temp + "\n");
                line = line + 1;
                bw.flush();
                bw.close();
            }   
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //现在是对每个亚种挑5个，输入文件是fa文件和model xml文件
    //输出的是特定的xml文件
    public void C39_get_xmlfile_forbeast(String infileS1,String infileS2,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp,temp2 = null;
            String key = null;
            String[] te = null;
            int line = 1,line2 = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);  
            BufferedReader br2 = IOUtils.getTextReader(infileS2);  
            String taxaNamelist[] = infileS1.toString().split("/");
            String filename = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
            HashMap<String, String> hashMapfa = new HashMap<String, String>();
            List<String> keylist = new ArrayList<String>();
             BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                if(line %2 == 1){
                    key = temp.substring(1, temp.length());
                }else{
                    hashMapfa.put(key, temp);
                    keylist.add(key);
                } 
                line = line +1;
            }  
            while((temp2 = br2.readLine())!=null){
                String newline = null;
                if(temp2.contains("A_2")){
                    temp2 = temp2.replaceAll("A_2", filename);
                }
                if(line2==5 || line2==6 || line2==7 || line2==8 || line2==9 || line2==10 || line2==11 || line2==12 || line2==13 || line2==14 || line2==15
                        || line2==16 || line2==17 || line2==18 || line2==19 || line2==20 || line2==21 || line2==22 || line2==23 || line2==24 || line2==25 
                        || line2==26 || line2==27 || line2==28 || line2==29 || line2==30){
                    for(int j = 0;j<keylist.size();j++){
                        if (temp2.contains(keylist.get(j))){
                            System.out.println(keylist.get(j));
                            String newvalue = "value=\"" + hashMapfa.get(keylist.get(j)) + "\"";
                            String old = "value=\"" + "[\\s\\S]*" + "\"";
                            newline = temp2.replaceAll(old, newvalue);
                            bw.write(newline);
                            bw.newLine();
                        }
                    }      
                }else{
                    bw.write(temp2);
                    bw.newLine();
                }  
                line2 = line2 +1;
            }                        
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //现在是对每个亚种挑1个,输入文件是fa文件和model xml文件
    //输出的是特定的xml文件
    public void C39_get_xmlfile_forbeast2(String infileS1,String infileS2,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp,temp2 = null;
            String key = null;
            String[] te = null;
            int line = 1,line2 = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);  
            BufferedReader br2 = IOUtils.getTextReader(infileS2);  
            String taxaNamelist[] = infileS1.toString().split("/");
            String filename = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
            HashMap<String, String> hashMapfa = new HashMap<String, String>();
            List<String> keylist = new ArrayList<String>();
             BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                if(line %2 == 1){
                    key = temp.substring(1, temp.length());
                }else{
                    hashMapfa.put(key, temp);
                    keylist.add(key);
                } 
                line = line +1;
            }  
            while((temp2 = br2.readLine())!=null){
                String newline = null;
                if(temp2.contains("species1line_6")){
                    temp2 = temp2.replaceAll("species1line_6", filename);
                }
                if(line2==5 || line2==6 || line2==7 || line2==8 || line2==9 || line2==10 || line2==11 || line2==12 || line2==13 || line2==14 || line2==15
                        || line2==16 || line2==17){
                    for(int j = 0;j<keylist.size();j++){
                        if (temp2.contains(keylist.get(j))){
                            System.out.println(keylist.get(j));
                            String newvalue = "value=\"" + hashMapfa.get(keylist.get(j)) + "\"";
                            String old = "value=\"" + "[\\s\\S]*" + "\"";
                            newline = temp2.replaceAll(old, newvalue);
                            bw.write(newline);
                            bw.newLine();
                        }
                    }      
                }else{
                    bw.write(temp2);
                    bw.newLine();
                }  
                line2 = line2 +1;
            }                        
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //现在是对每个亚种挑1个,输入文件是fa文件和model xml文件
    //输出的是特定的xml文件
    public void C39_get_xmlfile_forbeast3_A(String infileS1,String infileS2,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp,temp2 = null;
            String key = null;
            String[] te = null;
            int line = 1,line2 = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);  
            BufferedReader br2 = IOUtils.getTextReader(infileS2);  
            String taxaNamelist[] = infileS1.toString().split("/");
            String filename = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
            HashMap<String, String> hashMapfa = new HashMap<String, String>();
            List<String> keylist = new ArrayList<String>();
             BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                if(line %2 == 1){
                    key = temp.substring(1, temp.length());
                }else{
                    hashMapfa.put(key, temp);
                    keylist.add(key);
                } 
                line = line +1;
            }  
            while((temp2 = br2.readLine())!=null){
                String newline = null;
                if(temp2.contains("A_1line_4914")){
                    temp2 = temp2.replaceAll("A_1line_4914", filename);
                }
                if(line2==5 || line2==6 || line2==7 || line2==8 || line2==9 || line2==10 || line2==11 || line2==12 || line2==13 || line2==14 || line2==15
                        || line2==16 || line2==17 || line2==18 || line2==19 || line2==20 || line2==21 || line2==22 || line2==23 || line2==24 || line2==25 
                        || line2==26 ){
                    for(int j = 0;j<keylist.size();j++){
                        if (temp2.contains(keylist.get(j))){
                            System.out.println(keylist.get(j));
                            String newvalue = "value=\"" + hashMapfa.get(keylist.get(j)) + "\"";
                            String old = "value=\"" + "[\\s\\S]*" + "\"";
                            newline = temp2.replaceAll(old, newvalue);
                            bw.write(newline);
                            bw.newLine();
                        }
                    }      
                }else{
                    bw.write(temp2);
                    bw.newLine();
                }  
                line2 = line2 +1;
            }                        
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //现在是对每个亚种挑1个,输入文件是fa文件和model xml文件
    //输出的是特定的xml文件
    public void C39_get_xmlfile_forbeast3_B(String infileS1,String infileS2,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp,temp2 = null;
            String key = null;
            String[] te = null;
            int line = 1,line2 = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);  
            BufferedReader br2 = IOUtils.getTextReader(infileS2);  
            String taxaNamelist[] = infileS1.toString().split("/");
            String filename = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
            HashMap<String, String> hashMapfa = new HashMap<String, String>();
            List<String> keylist = new ArrayList<String>();
             BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                if(line %2 == 1){
                    key = temp.substring(1, temp.length());
                }else{
                    hashMapfa.put(key, temp);
                    keylist.add(key);
                } 
                line = line +1;
            }  
            while((temp2 = br2.readLine())!=null){
                String newline = null;
                if(temp2.contains("B_1line_4970")){
                    temp2 = temp2.replaceAll("B_1line_4970", filename);
                }
                if(line2==5 || line2==6 || line2==7 || line2==8 || line2==9 || line2==10 || line2==11 || line2==12 || line2==13 || line2==14 || line2==15
                        || line2==16 || line2==17 || line2==18 || line2==19 || line2==20 || line2==21 || line2==22 || line2==23 || line2==24 ){
                    for(int j = 0;j<keylist.size();j++){
                        if (temp2.contains(keylist.get(j))){
                            System.out.println(keylist.get(j));
                            String newvalue = "value=\"" + hashMapfa.get(keylist.get(j)) + "\"";
                            String old = "value=\"" + "[\\s\\S]*" + "\"";
                            newline = temp2.replaceAll(old, newvalue);
                            bw.write(newline);
                            bw.newLine();
                        }
                    }      
                }else{
                    bw.write(temp2);
                    bw.newLine();
                }  
                line2 = line2 +1;
            }                        
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //现在是对每个亚种挑1个,输入文件是fa文件和model xml文件
    //输出的是特定的xml文件
    public void C39_get_xmlfile_forbeast3_D(String infileS1,String infileS2,String outfileS){
        try {
            //BufferedWriter bw = null;
            String temp,temp2 = null;
            String key = null;
            String[] te = null;
            int line = 1,line2 = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);  
            BufferedReader br2 = IOUtils.getTextReader(infileS2);  
            String taxaNamelist[] = infileS1.toString().split("/");
            String filename = taxaNamelist[taxaNamelist.length-1].split("\\.")[0];
            HashMap<String, String> hashMapfa = new HashMap<String, String>();
            List<String> keylist = new ArrayList<String>();
             BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br.readLine())!=null){
                if(line %2 == 1){
                    key = temp.substring(1, temp.length());
                }else{
                    hashMapfa.put(key, temp);
                    keylist.add(key);
                } 
                line = line +1;
            }  
            while((temp2 = br2.readLine())!=null){
                String newline = null;
                if(temp2.contains("D_1line_4963")){
                    temp2 = temp2.replaceAll("D_1line_4963", filename);
                }
                if(line2==5 || line2==6 || line2==7 || line2==8 || line2==9 || line2==10 || line2==11 || line2==12 || line2==13 || line2==14 || line2==15
                        || line2==16 || line2==17 ){
                    for(int j = 0;j<keylist.size();j++){
                        if (temp2.contains(keylist.get(j))){
                            System.out.println(keylist.get(j));
                            String newvalue = "value=\"" + hashMapfa.get(keylist.get(j)) + "\"";
                            String old = "value=\"" + "[\\s\\S]*" + "\"";
                            newline = temp2.replaceAll(old, newvalue);
                            bw.write(newline);
                            bw.newLine();
                        }
                    }      
                }else{
                    bw.write(temp2);
                    bw.newLine();
                }  
                line2 = line2 +1;
            }                        
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    //这个基因的文件，42转成21
    public void FortheRealPos_gff(String infileS1,String infileS2,String outfileS){
        try{
            String temp = null;
            RowTable<String> genometable = new RowTable<>(infileS1);
            BufferedReader br = IOUtils.getTextReader(infileS2);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //bw.write(br.readLine() + "\n");
            while((temp = br.readLine()) != null){
                String[] tem = temp.split("\t");
                    int chr = Integer.valueOf(tem[0]);
                    String value1 = tem[2];
                    String value2 = tem[3];
                    String value3 = tem[4];
                    String pos = tem[1];
                    if(chr % 2 == 1){
                        String outchr = genometable.getCell(chr-1, 3);
                        //System.out.println(outchr);
                        bw.write(outchr + "\t" + pos + "\t" + value1 + "\t" + value2 + "\t" + value3 + "\n");
                    }else{
                        String outchr = genometable.getCell(chr-1, 3);
                        int outpos = Integer.valueOf(genometable.getCell(chr-1, 4)) + Integer.valueOf(pos);
                        bw.write(outchr + "\t" + outpos + "\t" + value1 + "\t" + value2 + "\t" + value3 +  "\n");
                    }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    //现在的方法是找到判断singlecopygene1_1.txt这个基因的长度是不是大于1000，要是不是大于1000，就不打印输出
//    String infileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/V1_1/GeneLulab1_1onlyGene.txt";
//    String infileS2 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1.txt";
//    String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/single_copy_gene/singlecopygene1_1_more1000.txt";
//    new split_time(infileS1,infileS2,outfileS);
    public void get_singlecopy_gene3_more1000(String infileS1,String infileS2, String outfileS){
        try {
            String temp = null;
            String temp1 = null;
            String[] te = null;
            int line = 1;
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedReader br2 = IOUtils.getTextReader(infileS2);
            BufferedWriter bw  = IOUtils.getTextWriter(outfileS);
            HashMap<String,Integer> hashMap1 = new HashMap<String,Integer>();
            while((temp = br.readLine()) != null){
                int length = Integer.valueOf(temp.split("\t")[3])-Integer.valueOf(temp.split("\t")[2]);
                hashMap1.put(temp.split("\t")[1],length);
            }
            while((temp1 = br2.readLine())!=null){
                String geneA = (temp1.split("\t")[2]).split("\\|")[1];
                String geneB = (temp1.split("\t")[3]).split("\\|")[1];
                String geneD = (temp1.split("\t")[4]).split("\\|")[1];
                if(hashMap1.get(geneA)>1000 && hashMap1.get(geneB)>1000 && hashMap1.get(geneD)>1000){
                    bw.write(temp1);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    //现在的方法是判断文件夹里面是不是有singlecopygene1_1_more1000.txt这个文件里面的基因，这是对文件夹的操作
    //必须要保证A，B，D三个基因组都是有的才可以
    //输入文件是singlecopygene1_1_more1000.txt这个文件和文件夹，输出的文件是三个文件夹，就是做cp
    //java -jar /data2/xuebo/Projects/Speciation/javaCode/C42_move1_1genefile.jar --file1 /data2/xuebo/Projects/Speciation/species_split_time/SwiftOrtho/singlecopygene1_1_more1000.txt --file2 /data1/home/xuebo/SRAssembler/assembly --out /data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/singlecopygene1_1_more1000_assembly.txt &
    public void C42_move1_1genefile(String infileS1,String infileS2,String outfileS){
        try{
            String temp = null;
            String temporder = null;
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            Set<String> genefile  = new HashSet<String>();
            File f = new File(infileS2);
            File[] tempList = f.listFiles();
            for (int i = 0; i < tempList.length; i++) {
                if (tempList[i].isFile()) {
                    //如果是文件，则不作处理
                }
                if (tempList[i].isDirectory()) {
                    String [] filename = tempList[i].toString().split("\\/");
                    genefile.add(filename[filename.length-1]);
                    System.out.println("fileName = " + filename[filename.length-1]);
                }
            }
            String dirA = "/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/A";
            String dirB = "/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/B";
            String dirD = "/data2/xuebo/Projects/Speciation/species_split_time/assemblygene/assembly_1_1/D";

//            String dirA = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/A";
//            String dirB = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/B";
//            String dirD = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/species_split_time/assembly_gene/D";
            while((temp = br.readLine()) != null){
                String geneA = (temp.split("\t")[2]).split("\\|")[1];
                String geneB = (temp.split("\t")[3]).split("\\|")[1];
                String geneD = (temp.split("\t")[4]).split("\\|")[1];
                if(!genefile.add(geneA) && !genefile.add(geneB) && !genefile.add(geneD)){ //表示组装好的基因里面是同时含有1_1基因的
                    bw.write(temp + "\n");
                    String  sourceDirA = infileS2.toString() + "/" + geneA;
                    copyDirectiory(sourceDirA, dirA + "/" + geneA);
                    String  sourceDirB = infileS2.toString() + "/" + geneB;
                    copyDirectiory(sourceDirB, dirB + "/" + geneB);
                    String  sourceDirD = infileS2.toString() + "/" + geneD;
                    copyDirectiory(sourceDirD, dirD + "/" + geneD);
                }else{
                    genefile.remove(geneA);
                    genefile.remove(geneB);
                    genefile.remove(geneD);
                }
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }


    public static void copyFile(File sourcefile, File targetFile) throws IOException {

        // 新建文件输入流并对它进行缓冲
        FileInputStream input = new FileInputStream(sourcefile);

        // 新建文件输出流并对它进行缓冲
        FileOutputStream out = new FileOutputStream(targetFile);
        BufferedOutputStream outbuff = new BufferedOutputStream(out);

        // 缓冲数组
        byte[] b = new byte[1024];
        int len = 0;
        while ((len = input.read(b)) != -1) {
            outbuff.write(b, 0, len);
        }

        //关闭文件
        outbuff.close();
        input.close();

    }

    public static void copyDirectiory(String sourceDir, String targetDir) throws IOException {

        // 新建目标目录
        (new File(targetDir)).mkdirs();

        // 获取源文件夹当下的文件或目录
        File[] file = (new File(sourceDir)).listFiles();

        for (int i = 0; i < file.length; i++) {
            if (file[i].isFile()) {
                // 源文件
                File sourceFile = file[i];
                // 目标文件
                File targetFile = new File(targetDir + File.separator + sourceFile.getName());
                copyFile(sourceFile, targetFile);

            }

            if (file[i].isDirectory()) {
                // 准备复制的源文件夹
                String dir1 = sourceDir + File.separator + file[i].getName();
                // 准备复制的目标文件夹
                String dir2 = targetDir + File.separator + file[i].getName();

                copyDirectiory(dir1, dir2);
            }
        }

    }








}
