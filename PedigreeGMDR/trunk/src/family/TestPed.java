package family;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;
import java.util.regex.Pattern;
import score.CalEngineException;
import edu.mit.wi.pedfile.PedFileException;
import family.pedigree.GMDRData;
import family.pedigree.GMDRPhenoFileException;
import family.pedigree.MDRPedFileException;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;

/**
 * It's a beta version of the PedGMDR, which produces nontransmitted genotypes and scores. TestPed is mainly developed
 * to realize the first step of PedGMDR, and the files produced by TestPed can be exported to GMDR.jar.
 * 
 * @author Guo-Bo Chen, Zhejiang University, chenguobo@gmail.com
 */
public class TestPed {

    /**
     * It reads parameters from the configuration file, which is given in the command line. When parameters are provided
     * correctly, the intermediate files will be produced by invoking those relevant classes. An IOException is thrown
     * when the number of parameters is not correct.
     * 
     * @param args
     *            a configuration file is required to provide parameters.
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
    	
		PedigreeParameter pr = new PedigreeParameter();
		if(args.length > 0) {
			pr.read(args[0]);
		} else {
			pr.setIsLouAlgorithm(true);
			pr.setUsingFounderGenotype(false);
			pr.setPedigreeFile("0.ped");
			pr.setPhenotypeFile("0.phe");
			pr.setConvertedPedigreeFile("Converted_0.ped");
			pr.setConvertedPhenotypeFile("Converted_0.phe");
			pr.setFamilyIDFile("Family_ID_0.txt");
			pr.setCovariateIndex("2");
			pr.setPhenotypeIndex("0");
			pr.setScoreBuildMethod(2);
			pr.setAdjustScore(true);
			pr.setScoreBuildWithFounder(true);
			pr.setReplication(10);
			pr.setSeed(10);
		}
		GMDRData.rnd = new Random(10);
		GMDRData GD = new GMDRData(pr.usingFounderGenotype, pr.isLouAlgorithm);
		AbstractGenoDistribution.rnd = new Random(pr.seed);
        boolean isPermutation = false; // 0 for all; 1 for affected only; 2 for unaffected only.

        GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());
        GD.RabinowitzApproach();
        GD.CreateTable(false);

		if (pr.getScoreBuildMethod() >= 0) {
			GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(), pr.getScoreBuildMethod(), pr.getScoreBuildWithFounder());
		} else {
			GD.fetchScore(pr.getPhenotypeIndex(), pr.getScoreBuildWithFounder());
		}
		GD.PrintGMDR(pr.getConvertedPedigreeFile(), pr.getConvertedPhenotypeFile(), isPermutation);
		isPermutation = true;
		for (int i = 0; i < pr.replication; i++) {
			String opfN = "Lou_" + Integer.toString(i) + ".txt";
			GD.PrintGMDR(opfN, null, isPermutation);
		}
    }
}
