package family;

import java.io.IOException;
import java.util.Random;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;

import score.CalEngineException;
import edu.mit.wi.pedfile.PedFileException;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.GMDRData;
import family.pedigree.GMDRPhenoFileException;
import family.pedigree.MDRPedFileException;
import mdr.data.DataFile;
import mdr.moore.AbstractMergeSearch;
import mdr.moore.LinearMergeSearch;
import mdr.GMDRParameter;

public class RunPedSimulation {

	public static void main (String[] args) throws IOException {

		PedigreeParameter pr = new PedigreeParameter();
		long seed = 10;
		if(args.length > 0) {
			pr.read(args[0]);
		} else {
			pr.setIsLouAlgorithm(false);
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
		boolean isRabinowitzProc = false;
		GMDRData.rnd = new Random(pr.seed);
		GMDRData GD = new GMDRData(pr.isLouAlgorithm, pr.usingFounderGenotype);
		AbstractGenoDistribution.rnd = new Random(pr.seed);
		GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());

		GD.RabinowitzApproach();
		GD.CreateTable(isRabinowitzProc);

		if (pr.getScoreBuildMethod() >= 0) {
			GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(), pr.getScoreBuildMethod(), 
					pr.getScoreBuildWithFounder());
		} else {
			GD.fetchScore(pr.getPhenotypeIndex(), pr.getScoreBuildWithFounder());
		}
		GD.PrintGMDR(pr.getConvertedPedigreeFile(), pr.getConvertedPhenotypeFile(), isRabinowitzProc);
		int[] scrIdx = {0};
		GMDRParameter gmdrPr = new GMDRParameter();
		gmdrPr.setInteractionEnd(2);
		gmdrPr.setInteractionFrom(2);
		gmdrPr.setInterval(5);
		gmdrPr.setMooreMDR(true);
		gmdrPr.setPartitionMethod(0);
		gmdrPr.setReplicationPermutation(10);
		gmdrPr.setSeed(seed);
		gmdrPr.setScoreIndex(scrIdx);
		DataFile mdrData = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD.getWorkingStatusTable(), GD.getTraitName(), GD.getWorkingScoreTable(), scrIdx);
        Subdivision sd = new Subdivision(gmdrPr.getInterval(), gmdrPr.getSeed(), mdrData);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(gmdrPr.getInterctionFrom(), gmdrPr.getInteractionEnd(), mdrData.getMarkerNum());
        cg.generateCombination();
        AbstractMergeSearch as;
        as = new LinearMergeSearch(mdrData, sd, cg, gmdrPr.getScoreIndex().length, mdrData.getOffset(), gmdrPr.isMooreMDR());
        for (int i = gmdrPr.getInterctionFrom(); i <= gmdrPr.getInteractionEnd(); i++) {
            as.search(i);
            as.summarise();
        }

        isRabinowitzProc = true;
		for (int i = 0; i < pr.replication; i++) {
			String opfN = "ParentsRabin_" + Integer.toString(i) + ".txt";
			GD.RabinowitzApproach();
			GD.CreateTable(isRabinowitzProc);
			GD.PrintGMDR(opfN, null, true);
		}
	}
}
