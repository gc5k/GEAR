package family;

import java.io.IOException;
import java.util.Random;
import java.util.ArrayList;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.GMDRData;
import mdr.data.DataFile;
import mdr.moore.AbstractMergeSearch;
import mdr.moore.LinearMergeSearch;
import mdr.GMDRParameter;

public class RunPedSimulation {

	public static void main(String[] args) throws IOException {

		for (int i = 0; i < 1; i++) {
			String PedFile = Integer.toString(i) + ".ped";
			String PhenoFile = Integer.toString(i) + ".phe";
			String CovertPedFile = "Converted_" + Integer.toString(i) + ".ped";
			String CovertPhenoFile = "Converted_" + Integer.toString(i)
					+ ".phe";
			String FamIDFile = "Family_ID_" + Integer.toString(i) + ".txt";
			PedigreeParameter pr = new PedigreeParameter();
			long seed = 10;
			if (args.length > 0) {
				pr.read(args[0]);
			} else {
				pr.setIsLouAlgorithm(false);
				pr.setUsingFounderGenotype(true);
				pr.setUsingChildrenGenotype(true);
				pr.setPedigreeFile(PedFile);
				pr.setPhenotypeFile(PhenoFile);
				pr.setConvertedPedigreeFile(CovertPedFile);
				pr.setConvertedPhenotypeFile(CovertPhenoFile);
				pr.setFamilyIDFile(FamIDFile);
				pr.setCovariateIndex("2");
				pr.setPhenotypeIndex("0");
				pr.setScoreBuildMethod(2);
				pr.setAdjustScore(true);
				pr.setScoreBuildWithFounder(true);
				pr.setScoreBuildWithChildren(true);
				pr.setReplication(10);
				pr.setSeed(10);
			}
			boolean isRabinowitzProc = false;
			GMDRData.rnd = new Random(pr.seed);
			GMDRData GD = new GMDRData(pr.usingFounderGenotype,
					pr.usingChildrenGenotype, pr.isLouAlgorithm);
			AbstractGenoDistribution.rnd = new Random(pr.seed);
			GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());

//			GD.RabinowitzApproach();
			GD.CreateTable(isRabinowitzProc);

			if (pr.getScoreBuildMethod() >= 0) {
				GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(),
						pr.isAdjustScore(), pr.getScoreBuildMethod(), pr
								.getScoreBuildWithFounder(), pr
								.getScoreBuildWithChildren());
			} else {
				GD.fetchScore(pr.getPhenotypeIndex(), pr
						.getScoreBuildWithFounder());
			}
			GD.CreateWorkingTable(isRabinowitzProc);
			// GD.PrintGMDR(pr.getConvertedPedigreeFile(),
			// pr.getConvertedPhenotypeFile(), isRabinowitzProc);
			int[] scrIdx = { 0 };
			GMDRParameter gmdrPr = new GMDRParameter();
			gmdrPr.setInteractionEnd(2);
			gmdrPr.setInteractionFrom(2);
			gmdrPr.setInterval(5);
			gmdrPr.setMooreMDR(true);
			gmdrPr.setPartitionMethod(0);
			gmdrPr.setReplicationPermutation(10);
			gmdrPr.setSeed(seed);
			gmdrPr.setScoreIndex(scrIdx);
			DataFile mdrData = new DataFile(GD.getMarkerName(), GD
					.getWorkingGenoTable(), GD.getWorkingStatusTable(), GD
					.getTraitName(), GD.getWorkingScoreTable(), scrIdx);
			Subdivision sd = new Subdivision(gmdrPr.getInterval(), gmdrPr
					.getSeed(), mdrData);
			sd.RandomPartition();

			CombinationGenerator cg = new CombinationGenerator(gmdrPr
					.getInterctionFrom(), gmdrPr.getInteractionEnd(), mdrData
					.getMarkerNum());
			cg.generateCombination();
			AbstractMergeSearch as;
			as = new LinearMergeSearch(mdrData, sd, cg,
					gmdrPr.getScoreIndex().length, mdrData.getOffset(), gmdrPr
							.isMooreMDR());
			for (int j = gmdrPr.getInterctionFrom(); j <= gmdrPr.getInteractionEnd(); j++) {
				as.search(j);
				as.summarise();
				int[] bm = as.getBestModel(j, scrIdx[0]);
				GD.SetChosenMarker(bm);
				isRabinowitzProc = false;
				System.out.println("permutation test");
				for (int k = 0; k < pr.replication; k++) {
					GD.RabinowitzApproach();
					GD.CreateTable(isRabinowitzProc);
					String opfN = "ParentsRabin_" + Integer.toString(k) + ".txt";
					GD.CreateWorkingTable(isRabinowitzProc);
					DataFile mdrData1 = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD.getWorkingStatusTable(),
							GD.getTraitName(), GD.getWorkingScoreTable(), scrIdx);

					CombinationGenerator cg1 = new CombinationGenerator(j, j, mdrData1
							.getMarkerNum());
					cg1.generateCombination();
					AbstractMergeSearch as1;
					Subdivision sd1 = new Subdivision(gmdrPr.getInterval(), gmdrPr
							.getSeed(), mdrData1);
					sd1.RandomPartition();
					as1 = new LinearMergeSearch(mdrData1, sd1, cg1,
							gmdrPr.getScoreIndex().length, mdrData1.getOffset(), gmdrPr
									.isMooreMDR());
					as1.search(j);
					as1.summarise();
				}
			}
		}
	}
}
