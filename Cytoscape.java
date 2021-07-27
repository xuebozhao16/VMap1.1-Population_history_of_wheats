/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

//import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class Cytoscape {
    public Cytoscape(String infileS,String outfileS){
        this.getStandardphy_file(infileS, outfileS);
    }
    
    //这个方法是为了产生标准的phy文件
    public void getStandardphy_file(String infileS,String outfileS){
        try{
            String temp = null;
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String line1 = br.readLine();
            bw.write(line1 + "\n");
            while((temp = br.readLine()) != null){
                String[] tem = temp.split(" ");
                int length = 10 - tem[0].length();
                bw.write(tem[0] + multipleSpaces(length) + tem[1] + "\n");               
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    public String multipleSpaces(int n){
        String output = ""; 
        for(int i=0; i<n; i++)
           output += " ";
        return output;
    }
}
