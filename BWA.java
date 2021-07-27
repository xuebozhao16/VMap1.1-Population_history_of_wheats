/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import pgl.infra.utils.IOUtils;
//import utils.IOUtils;

import static pgl.infra.utils.IOUtils.listFilesEndsWith;
import static pgl.infra.utils.IOUtils.listRecursiveFiles;
//import static utils.IOUtils.getTextReader;
//import static utils.IOUtils.getTextWriter;
//import static utils.IOUtils.listFilesEndsWith;
//import static utils.IOUtils.listRecursiveFiles;

/**
 *
 * @author xuebozhao
 */
public class BWA {
    public BWA(String infileS1,String infileS2,String outfileS){
        this.bwa_file(infileS1, infileS2, outfileS);
    }
    public BWA(String infileS1,String infileS2,String outfileS1,String outfileS2){
        this.C1_testBWA(infileS1, infileS2, outfileS1, outfileS2);
    }
    
    //这个方法是生成run bwa所需要的文件
    // GetWheatCDSsequence  new BWA(infileS1,infileS2,outfileS);
//        String infileS1 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/bwamap/bwa_file1.txt";
//        String infileS2 = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/bwamap/test/AT19767A";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/bwamap/test/bwatest.txt";
//        new BWA(infileS1,infileS2,outfileS);
    public void bwa_file(String infileS1,String infileS2,String outfileS){
        try{    
            String temp = null;
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            HashMap<String, String> hashMap1 = new HashMap<String, String>();
            while((temp = br.readLine()) != null){
                    String[] tem = temp.split("\t");
                    String key = tem[0] + "A";
                    String value = tem[1] + "__" + tem[2];
                    hashMap1.put(key, value);
                    
            }
            //BufferedReader NGSFile = null;
            File f = new File(infileS2);
            String infileS2name = infileS2.substring(infileS2.length()-8, infileS2.length());
            System.out.println(infileS2name + "\n");
            File[] fs = listRecursiveFiles(f);
            File[] sub = listFilesEndsWith(fs, ".clean.fq.gz");
            for(File fi:sub){  
                String out = fi.toString();
                if (out.endsWith("_1.clean.fq.gz")){
                    System.out.println(out + "\n");
                    bw.write(out + "\t");
                }
                if (out.endsWith("_2.clean.fq.gz")){
                    bw.write(out + "\t" + hashMap1.get(infileS2name).split("__")[0] + "\t" 
                        + hashMap1.get(infileS2name).split("__")[1] + "\n");
                }    

            }             
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }


    ///data3/wgs/fastq/xuebo_S20/AT19767A/190923_I013_V300031142_L3_AE07750767-510 随机挑100K的reads 要成对哦
    //这个和RandomLine.java这个代码不一样，因为这个要记录文件的位置，下一个文件也要挑这些
    //这个文件122085404行，，122085404/4 = 30521351 ,,在30521351(0-30521350)里面随机抽100,000,得到的位点*4+1即可
    public void C1_testBWA(String infileS1,String infileS2,String outfileS1,String outfileS2){
        try{    
            String temp = null;
            String temp2 = null;
            BufferedReader br1 = null; 
            BufferedReader br2 = null; 
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
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
            Set SetNum = new HashSet();
            //产生随机数
            Random rand = new Random();
            for(int i=0; i<5; i++) {
                //产生的随机数放到SetNum里面
                System.out.println(rand.nextInt(10) + 1);
                int AA = rand.nextInt(10) + 1;
                SetNum.add(AA*4+1);
            }
            int line1 = 1;
            while((temp = br1.readLine()) != null){
                    if(!SetNum.add(line1)){
                        String name = temp.replace("/1", "/2");
                        Set1.add(name);
                        bw1.write(temp + "\n");
                        bw1.write(br1.readLine() + "\n");
                        bw1.write(br1.readLine() + "\n");
                        bw1.write(br1.readLine() + "\n");
                        line1 = line1 + 4;
                    }else{
                       SetNum.remove(line1);
                    }   
                line1 = line1 + 4;
            }                     
            bw1.flush();
            bw1.close();
            while((temp2 = br2.readLine()) != null){
                if(!Set1.add(temp2)){
                    bw2.write(temp2 + "\n");
                    bw2.write (br2.readLine() + "\n");
                    bw2.write (br2.readLine() + "\n");
                    bw2.write (br2.readLine() + "\n");
                }else{
                    Set1.remove(temp2);
                }
            }
            bw2.flush();
            bw2.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
