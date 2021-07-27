package speciation;

import com.sun.org.apache.xerces.internal.xs.StringList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;

public class IBS_distance {
    public IBS_distance(String infileS1,String infileS2,String outfileS){
        this.C40_getIBS_distance(infileS1, infileS2, outfileS);
    }
    public IBS_distance(String infileS1,String outfileS){
        this.C41_getIBS_distance2(infileS1, outfileS);
    }


    //这个function的目的是对两个vcf进行的处理
    public void C40_getIBS_distance(String infileS1,String infileS2,String outfileS){
        try{
            GenotypeGrid g1 = new GenotypeGrid(infileS1, GenoIOFormat.VCF);
            GenotypeGrid g2 = new GenotypeGrid(infileS2, GenoIOFormat.VCF);
            //GenotypeGrid f=new GenotypeGrid("",GenoIOFormat.Binary);
            GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1,g2);
            SumTaxaDivergence std = new SumTaxaDivergence(g);
            std.writeDxyMatrix(outfileS, IOFileFormat.Text);
            //g.getIBSDistanceMatrix();
            //System.out.println(g.getIBSDistanceMatrix());
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    
    //这个function的目的是对一个vcf进行的处理
    public void C41_getIBS_distance2(String infileS1,String outfileS){
        try{
            GenotypeGrid g1 = new GenotypeGrid(infileS1, GenoIOFormat.VCF);
            SumTaxaDivergence std = new SumTaxaDivergence(g1);
            std.writeDxyMatrix(outfileS, IOFileFormat.Text);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
