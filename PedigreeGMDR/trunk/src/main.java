
import java.io.*;
import java.util.ArrayList;
import java.util.regex.Pattern;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;

import mdr.moore.LinearMergeSearch;
import mdr.moore.Permutation;
import mdr.moore.PermutationResult;
import mdr.SavedModels;
import mdr.moore.AbstractMergeSearch;
import mdr.data.DataFile;
import mdr.DataFileException;
import mdr.GMDRParameter;

/**
 *
 * @author Guo-Bo Chen
 */
public class main {

    public static void main(String[] args) throws IOException {
    	GMDRParameter pr = new GMDRParameter();
    	if (args.length > 0) {
    		pr.read(args[0]);
    	} else {
    		pr.setMarkerFile("ped0.txt");
    		pr.setPhenotypeFile("phe0.txt");
    		pr.setInteractionFrom(1);
    		pr.setInteractionEnd(2);
    		pr.setScoreIndex("0");
    		pr.setInterval(10);    		
    		pr.setSeed(90);
    		pr.setPartitionMethod(0);
    		pr.setMooreMDR(true);
    		pr.setReplicationPermutation(100);
    		pr.setSearchMethod(0);

    	}

//        pFile = null;
        DataFile dr = new DataFile(pr.getMarkerFile(), pr.getPhenotypeFile(), pr.getScoreIndex());

        Subdivision sd = new Subdivision(pr.getInterval(), pr.getSeed(), dr);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(pr.getInterctionFrom(), pr.getInteractionEnd(), dr.getMarkerNum());
        cg.generateCombination();
        AbstractMergeSearch as = new LinearMergeSearch(dr, sd, cg, pr.getScoreIndex().length, dr.getOffset(), pr.isMooreMDR());
        for (int i = pr.getInterctionFrom(); i <= pr.getInteractionEnd(); i++) {
            as.search(i);
            as.summarise();
        }

        for (int i = pr.getInterctionFrom(); i <= pr.getInteractionEnd(); i++) {
            SavedModels sm = as.getBestSavedModelAtOrder(new Integer(i));
            Permutation pm = new Permutation(dr, sm, dr.getOffset(), sd, pr.getReplicationPermutation());
            pm.evaluate();
            PermutationResult permu_res = pm.getResult(0);
            permu_res.summarise();
            System.out.println(permu_res);
        }
    }
}
