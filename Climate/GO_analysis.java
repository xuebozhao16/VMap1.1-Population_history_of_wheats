/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
//import utils.IOUtils;
import xuebo.analysis.annotation.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class GO_analysis {
    public GO_analysis(String infileS1,String infileS2, String outfileS){
        this.forGOandGENEforspecifiedGene(infileS1, infileS2, outfileS);
    }
    

    //这个方法是进行GO分析，得到文件1：基因号和GO对应的文件 文件2：GO和GO的注释信息
    public void forGOandGENEforspecifiedGene(String infileS1,String infileS2, String outfileS){
        try{
            String tempgene = null;
            String temp = null;
            BufferedReader brspecifiedGene = IOUtils.getTextReader(infileS1);
            BufferedReader brGoFile = IOUtils.getTextReader(infileS2);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS);
            //BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            //BufferedWriter bw3 = IOUtils.getTextWriter(outfileS3);
            Set specifiedGene = new HashSet();
            Set GoAnnotation1 = new HashSet();
            Set GoAnnotation11 = new HashSet();
            Set GoAnnotation2 = new HashSet();
            while((tempgene = brspecifiedGene.readLine()) != null){
                specifiedGene.add(tempgene);               
            }
            brGoFile.readLine();
            while((temp = brGoFile.readLine()) !=null){
                String[] tem = temp.split("\t");
                if(!specifiedGene.add(tem[0])){
                    GoAnnotation1.add(tem[2] + "_" + tem[0]);
                    //GoAnnotation1.add(tem[4] + "_" + tem[0]);
                    GoAnnotation11.add(tem[0]);
                    GoAnnotation2.add(tem[2]+"_" +tem[3]);
                    //GoAnnotation2.add(tem[4]+"_" +tem[5]);
                }else{
                    specifiedGene.remove(tem[0]);
                }
            }
            for (Object str1 : GoAnnotation1) {
                String GOgene = str1.toString();
                //bw1.write(GOgene.split("_")[0] + "," + GOgene.split("_")[1] + "\n");
                bw1.write(GOgene.split("_")[1] + "\t" + GOgene.split("_")[0] + "\n");
            }           
            bw1.flush();
            //bw2.flush();
            //bw3.flush();
            bw1.close();
            //bw2.close();
            //bw3.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
