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
		Report report = new Report();
		int replication = args.length > 0 ? Integer.parseInt(args[0]) : 5; 
		for (int i = 0; i < replication; i++) {
			String PedFile = Integer.toString(i) + ".ped";
			String PhenoFile = Integer.toString(i) + ".phe";
			String CovertPedFile = "Converted_" + Integer.toString(i) + ".ped";
			String CovertPhenoFile = "Converted_" + Integer.toString(i)
					+ ".phe";
			String FamIDFile = "Family_ID_" + Integer.toString(i) + ".txt";
			PedigreeParameter pr = new PedigreeParameter();
			long seed = 10;
			if (args.length > 1) {
				pr.read(args[1]);
				pr.setPedigreeFile(PedFile);
				pr.setPhenotypeFile(PhenoFile);
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
			if (args.length > 2) {
				gmdrPr.read(args[2]);
			} else {
				gmdrPr.setInteractionEnd(2);
				gmdrPr.setInteractionFrom(2);
				int[] scrIdx = { 0 };
				gmdrPr.setScoreIndex(scrIdx);
				gmdrPr.setInterval(5);
				gmdrPr.setSeed(seed);
				gmdrPr.setPartitionMethod(0);
				gmdrPr.setMooreMDR(true);
				gmdrPr.setReplicationPermutation(10);
				gmdrPr.setSearchMethod(0);
			}

			boolean isRabinowitzProc = false;
			GMDRData.rnd = new Random(pr.seed + i);
			GMDRData GD = new GMDRData(pr.usingFounderGenotype,
					pr.usingChildrenGenotype, pr.isLouAlgorithm);
			AbstractGenoDistribution.rnd = new Random(pr.seed + 1);
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
					.getSeed() + i, mdrData);
			sd.RandomPartition();

			CombinationGenerator cg = new CombinationGenerator(gmdrPr
					.getInterctionFrom(), gmdrPr.getInteractionEnd(), mdrData
					.getMarkerNum());
			cg.generateCombination();
			LinearMergeSearch as;
			as = new LinearMergeSearch(mdrData, sd, cg,
					gmdrPr.getScoreIndex().length, mdrData.getOffset(), gmdrPr
							.isMooreMDR());
			for (int j = gmdrPr.getInterctionFrom(); j <= gmdrPr.getInteractionEnd(); j++) {
				as.search(j);
				report.NewRound(as.getBestModelKey(j, 0));
				double[][] stats = as.singleBest(as.getBestModelKey(j, 0));
				report.Add_test_statistic(stats[0]);
				isRabinowitzProc = false;
				for (int k = 0; k < pr.replication; k++) {
					isRabinowitzProc = true;
					GD.RabinowitzApproach();
					GD.CreateTable(isRabinowitzProc);
					GD.CreateWorkingTable(isRabinowitzProc);
					DataFile mdrData1 = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD.getWorkingStatusTable(),
							GD.getTraitName(), GD.getWorkingScoreTable(), gmdrPr.getScoreIndex());
					CombinationGenerator cg1 = new CombinationGenerator(j, j, mdrData1
							.getMarkerNum());
					cg1.generateCombination();
					LinearMergeSearch as1;
					Subdivision sd1 = new Subdivision(gmdrPr.getInterval(), gmdrPr
							.getSeed()+ i * pr.replication + k, mdrData1);
					sd1.RandomPartition();
					as1 = new LinearMergeSearch(mdrData1, sd1, cg1,
							gmdrPr.getScoreIndex().length, mdrData1.getOffset(), gmdrPr
									.isMooreMDR());
					as1.search(j);
					double[][] stats1 = as1.singleBest(as1.getBestModelKey(j, 0));
					report.Add_null_test_statistic(stats1[0]);
				}
				report.RoundSummary();
			}
		}
		System.out.println(report);
	}
}
