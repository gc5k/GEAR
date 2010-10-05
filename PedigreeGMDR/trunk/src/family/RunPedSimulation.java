package family;

import java.io.IOException;
import java.util.Random;

import score.CalEngineException;
import edu.mit.wi.pedfile.PedFileException;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.GMDRData;
import family.pedigree.GMDRPhenoFileException;
import family.pedigree.MDRPedFileException;

public class RunPedSimulation {

	public static void main (String[] args) throws IOException {

		PedigreeParameter pr = new PedigreeParameter();
		if(args.length > 0) {
			pr.read(args[0]);
		} else {
			pr.setPedigreeFile("0.ped");
			pr.setPhenotypeFile("0.phe");
			pr.setConvertedPedigreeFile("Converted_0.ped");
			pr.setConvertedPhenotypeFile("Converted_0.phe");
			pr.setFamilyIDFile("Family_ID_0.txt");
			pr.setCovariateIndex("1");
			pr.setPhenotypeIndex("0");
			pr.setScoreBuildMethod(1);
			pr.setAdjustScore(true);
			pr.setScoreBuildWithFounder(true);
			pr.setReplication(10);
		}

		GMDRData GD = new GMDRData(true);

		try {
			GD.InitialPedFile(pr.getPedigreeFile());// initial Pedfile
			GD.InitialPhenoFile(pr.getPhenotypeFile());// initial phenotype
		} catch (MDRPedFileException E) {
			E.printStackTrace(System.err);
		} catch (GMDRPhenoFileException e) {
			System.err.println("Phenotype File Exception.");
		}
		AbstractGenoDistribution.rnd = new Random(2);
		GD.Match();
		GD.RabinowitzApproach();
		GD.realCreateTableWithParents();
		try {
			if (pr.getScoreBuildMethod() >= 0) {
				GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(), pr.getScoreBuildMethod(), 
						pr.getScoreBuildWithFounder());
			} else {
				GD.fetchScore(pr.phe_idx[0]);
			}
			GD.RabinowitzPrintGMDRwithParents(pr.converted_ped_file, pr.converted_phe_file, false);
		} catch (CalEngineException e) {
			e.printStackTrace(System.err);
		} catch (Exception e) {
			e.printStackTrace(System.err);
		}

		
		
		for (int i = 0; i < pr.replication; i++) {
			String opfN = "ParentsRabin_" + Integer.toString(i) + ".txt";
			GD.RabinowitzApproach();
			GD.RabinowitzCreateTableWithParents();
			try {
				GD.RabinowitzPrintGMDRwithParents(opfN, null, true);
			} catch (CalEngineException e) {
				e.printStackTrace(System.err);
			} catch (Exception e) {
				e.printStackTrace(System.err);
			}
		}
	}
}
