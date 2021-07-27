/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package speciation;

import gnu.trove.list.array.TDoubleArrayList;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import pgl.infra.dna.FastaByte;
import pgl.infra.dna.FastaRecordByte;

import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class vcf_QualityControl {
//    public vcf_QualityControl(String infileS, String outfileS1,String outfileS2){
//        this.C36_checkQuality(infileS, outfileS1, outfileS2);
//    }
    public vcf_QualityControl(String infileS1, String infileS2,String outfileS){
            this.C37_getReference_seq(infileS1,infileS2,outfileS);
    }
    public void C36_checkQuality (String infileS, String outfileS1,String outfileS2){
        try {
            //现在是开始计算
            TDoubleArrayList missingSite = new TDoubleArrayList(); // calculation 1
            TDoubleArrayList hetSite = new TDoubleArrayList(); // calculation 2
            TDoubleArrayList maf = new TDoubleArrayList(); // calculation 3
            TDoubleArrayList missingTaxon = new TDoubleArrayList(); // calculation 4
            TDoubleArrayList hetTaxon = new TDoubleArrayList(); // calculation 5      
            GenotypeGrid gt = new GenotypeGrid(infileS,GenoIOFormat.VCF_GZ);
            for (int i = 0; i < gt.getSiteNumber(); i=i+1) {
                //System.out.println(gt.getChromosome(i));
                //System.out.println(gt.getPosition(i));
                missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
                hetSite.add(gt.getHeterozygousProportionBySite(i));
                maf.add(gt.getMinorAlleleFrequency(i));
            }
            for (int i = 0; i < gt.getTaxaNumber(); i++) {
                missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
                hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
            }
            //现在是开始输出bw1
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            bw1.write("Chr\tPos\tHeterozygousProportion\tMissingRate\tMaf");
            bw1.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < missingSite.size(); i++) {
                sb.setLength(0);
                sb.append(gt.getChromosome(i)).append("\t").append(gt.getPosition(i)).append("\t").
                        append(hetSite.get(i)).append("\t").append(missingSite.get(i)).append("\t").append(maf.get(i));
                bw1.write(sb.toString());
                bw1.newLine();
            }
            //现在是开始输出bw2
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            bw2.write("Taxa\tHeterozygousProportion\tMissRate");
            bw2.newLine();
            StringBuilder sb2 = new StringBuilder();
            for (int i = 0; i < gt.getTaxaNumber(); i++) {
                sb2.setLength(0);
                sb2.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i)).append("\t").append(missingTaxon.get(i));
                bw2.write(sb2.toString());
                bw2.newLine();
            }
            bw1.flush();
            bw2.flush();
            bw1.close();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    //现在这个方法是3 bit的操作，截取基因组上面的指定位置信息的fasta
    //输入文件是小麦的基因组，bed文件，得到的是序列的信息
    public void C37_getReference_seq (String infileS1, String infileS2,String outfileS){
        try {
            String temp2 = null;
            String chrseq  = null;
            BufferedReader br2 = IOUtils.getTextReader(infileS2);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef");
            bw.newLine();
            //Set vcfset = new HashSet();
            FastaByte chrFa=new FastaByte(infileS1);
            while((temp2 = br2.readLine())!= null){
                int chr = Integer.valueOf(temp2.split("\t")[0]);
                int vcfpos = Integer.valueOf(temp2.split("\t")[1]);
                int chrIndex = chrFa.getIndexByName(String.valueOf(chr));
                chrseq = chrFa.getSeq(chrIndex);
                //System.out.println(chrseq);
                String ref = chrseq.substring(vcfpos-1,vcfpos);
                //gt.getSequence();
                //System.out.println(ref);
                bw.write(chr + "\t" + vcfpos + "\t" + ref + "\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }




    
    
    
//    private void Aoyue_checkQuality (String vcfDirS, String outDirS, String genomeType) {
//        File outDir = new File (outDirS);
//        outDir.mkdir();
//        int fN = 2;
//        int size = 20000; //抽样数量
//        List<File> fList = IOUtils.getFileListInDir(vcfDirS);
//        Collections.sort(fList);
//
//        GenotypeGrid[] gts = new GenotypeGrid[fN];
//        int totalSiteCount = 0;
//        for (int i = 0; i < gts.length; i++) { //对genotypeTable 进行初始化
//            gts[i] = new GenotypeGrid(fList.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
//            totalSiteCount+=gts[i].getSiteNumber(); //获取 fN 个文件的总位点数
//        }
//        GenotypeOperation.mergeGenotypesBySite(gts[0], gts[1]); // 合并两个VCF文件
//        GenotypeGrid gt = gts[0];
//        TDoubleArrayList missingSite = new TDoubleArrayList(); // calculation 1
//        TDoubleArrayList hetSite = new TDoubleArrayList(); // calculation 2
//        TDoubleArrayList maf = new TDoubleArrayList(); // calculation 3
//        TDoubleArrayList missingTaxon = new TDoubleArrayList(); // calculation 4
//        TDoubleArrayList hetTaxon = new TDoubleArrayList(); // calculation 5
//
//        int step = gt.getSiteNumber()/size;
//        for (int i = 0; i < gt.getSiteNumber(); i+=step) {
//            missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
//            hetSite.add(gt.getHeterozygousProportionBySite(i));
//            maf.add(gt.getMinorAlleleFrequency(i));
//        }
//        for (int i = 0; i < gt.getTaxaNumber(); i++) {
//            missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
//            hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
//        }
//
//        String siteQCfileS = new File(outDirS,genomeType + "_site_QC.txt.gz").getAbsolutePath();
//        String taxaQCFileS = new File (outDirS, genomeType+"_taxa_QC.txt.gz").getAbsolutePath();
//
//        try {
//            BufferedWriter bw = IOUtils.writeFile(siteQCfileS);
//            bw.write("GenomeType\tHeterozygousProportion\tMissingRate\tMaf");
//            bw.newLine();
//            StringBuilder sb = new StringBuilder();
//            for (int i = 0; i < missingSite.size(); i++) {
//                sb.setLength(0);
//                sb.append(genomeType).append("\t").append(hetSite.get(i)).append("\t").append(missingSite.get(i)).append("\t").append(maf.get(i));
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//
//
//        try {
//            BufferedWriter bw = AoFile.writeFile(taxaQCFileS);
//            bw.write("Taxa\tHeterozygousProportion\tMissRate\tGenomeType");
//            bw.newLine();
//            StringBuilder sb = new StringBuilder();
//            for (int i = 0; i < gt.getTaxaNumber(); i++) {
//                sb.setLength(0);
//                sb.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i)).append("\t").append(missingTaxon.get(i)).append("\t").append(genomeType);
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
}
