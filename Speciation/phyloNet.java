/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

/**
 *
 * @author xuebozhao
 */
//这个文件是产生phyloNet相关的代码，包括怎样产生数据和结果文件的分析展示
public class phyloNet {
    public phyloNet(String infileS,String outfileS){
        this.getphyloNet_tree(infileS, outfileS);
    }
    
    public void getphyloNet_tree(String infileS,String outfileS){
        try{
            String temp = null;
            int i = 0;
            BufferedReader br1 = IOUtils.getTextReader(infileS); 
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while((temp = br1.readLine())!= null){
                //String tem[] = temp.split("\t");
                bw.write("Tree gt" + i + "=" + temp + "\n");
                i++;
            }

                bw.flush();
                bw.close();               
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
