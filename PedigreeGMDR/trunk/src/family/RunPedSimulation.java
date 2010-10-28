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
import family.report.Report;

public class RunPedSimulation {

	public static void main(String[] args) throws IOException {
		String truemodel = "0,3";
		Report report = new Report();
		for (int i = 0; i < 5; i++) {

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
				pr.setUsingChildrenGenotype(false);
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
				pr.setScoreBuildWithChildren(false);
				pr.setReplication(10);
				pr.setSeed(10);
			}
			GMDRParameter gmdrPr = new GMDRParameter();
			if (args.length > 1) {
				gmdrPr.read(args[1]);
			} else {
				int[] scrIdx = { 0 };
				gmdrPr.setInteractionEnd(2);
				gmdrPr.setInteractionFrom(2);
				gmdrPr.setInterval(5);
				gmdrPr.setMooreMDR(true);
				gmdrPr.setPartitionMethod(0);
				gmdrPr.setReplicationPermutation(10);
				gmdrPr.setSeed(seed);
				gmdrPr.setScoreIndex(scrIdx);				
			}

			boolean isRabinowitzProc = false;
			GMDRData.rnd = new Random(pr.seed);
			GMDRData GD = new GMDRData(pr.usingFounderGenotype,
					pr.usingChildrenGenotype, pr.isLouAlgorithm);
			AbstractGenoDistribution.rnd = new Random(pr.seed);
			GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());

			GD.CreateTable(isRabinowitzProc);
			if (pr.getScoreBuildMethod() >= 0) {
				GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(),
						pr.isAdjustScore(), pr.getScoreBuildMethod(),
						pr.getScoreBuildWithFounder(),
						pr.getScoreBuildWithChildren());
			} else {
				GD.fetchScore(pr.getPhenotypeIndex(), pr
						.getScoreBuildWithFounder());
			}
			GD.CreateWorkingTable(isRabinowitzProc);

			DataFile mdrData = new DataFile(GD.getMarkerName(), GD
					.getWorkingGenoTable(), GD.getWorkingStatusTable(), GD
					.getTraitName(), GD.getWorkingScoreTable(), gmdrPr.getScoreIndex());
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
				report.NewRound(as.getBestModelKey(j, 0));
				int[] bm = as.getBestModel(j, 0);
				GD.SetChosenMarker(bm);
				isRabinowitzProc = false;
				for (int k = 0; k <= pr.replication; k++) {
					GD.RabinowitzApproach();
					GD.CreateTable(isRabinowitzProc);
					GD.CreateWorkingTable(isRabinowitzProc);
					DataFile mdrData1 = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD.getWorkingStatusTable(),
							GD.getTraitName(), GD.getWorkingScoreTable(), gmdrPr.getScoreIndex());

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
					if (k == 0) {
						double[] s = as1.getStats(0, j);
						report.Add_test_statistic(s);
					} else {
						double[] s = as1.getStats(0, j);
						report.Add_null_test_statistic(s);
					}
					isRabinowitzProc = true;
				}
				report.RoundSummary();		
			}
		}
		System.out.println(report);
	}
}
