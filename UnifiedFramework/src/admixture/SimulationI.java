package admixture;

import java.io.IOException;
import java.util.Random;

import mdr.GMDRParameter;
import mdr.algorithm.CombinationGenerator;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import mdr.moore.AbstractMergeSearch;
import mdr.moore.HeteroLinearMergeSearch;
import mdr.moore.LinearMergeSearch;
import family.PedigreeParameter;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.ChenAlgorithm;

public class SimulationI {

	public static void main(String[] args) throws IOException {
		int replication = args.length > 0 ? Integer.parseInt(args[0]) : 1;
		for (int i = 0; i < replication; i++) {
			String PedFile = Integer.toString(i) + "L_ped.txt";
			String PhenoFile = Integer.toString(i) + "L_phe.txt";
			String CovertPedFile = "Converted_" + Integer.toString(i) + ".ped";
			String CovertPhenoFile = "Converted_" + Integer.toString(i) + ".phe";
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
				pr.setUsingChildrenGenotype(true);
				pr.setPedigreeFile(PedFile);
				pr.setPhenotypeFile(PhenoFile);
				pr.setConvertedPedigreeFile(CovertPedFile);
				pr.setConvertedPhenotypeFile(CovertPhenoFile);
				pr.setFamilyIDFile(FamIDFile);
				pr.setCovariateIndex("2");
				pr.setPhenotypeIndex("1");
				pr.setScoreBuildMethod(0);
				pr.setAdjustScore(false);
				pr.setScoreBuildWithFounder(true);
				pr.setScoreBuildWithChildren(true);
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
				gmdrPr.setReplicationPermutation(5);
				gmdrPr.setSearchMethod(0);
			}

			ChenAlgorithm.rnd = new Random(pr.getSeed() + i);
			ChenAlgorithm GD = new ChenAlgorithm();
			AbstractGenoDistribution.rnd = new Random(pr.getSeed() + 1);
			GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());

			if (pr.getScoreBuildMethod() >= 0) {
				GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(), pr.getScoreBuildMethod(), pr
						.getScoreBuildWithFounder(), pr.getScoreBuildWithChildren());
			} else {
				GD.fetchScore(pr.getPhenotypeIndex(), pr.getScoreBuildWithFounder());
			}

			DataFile mdrData = new DataFile(GD.getMarkerName(), GD.getGenotype(), 
					GD.getStatue(), GD.getScoreName(), GD.getScore2());
			Subdivision sd = new Subdivision(gmdrPr.getInterval(), gmdrPr.getSeed() + i, mdrData.size());
			sd.RandomPartition();

			AbstractMergeSearch as = new HeteroLinearMergeSearch(mdrData, sd, gmdrPr.isMooreMDR());
			for (int j = gmdrPr.getInterctionFrom(); j <= gmdrPr.getInteractionEnd(); j++) {
				CombinationGenerator cg = new CombinationGenerator(j,
						j,	mdrData.getNumMarker());
				cg.generateCombination();
				as.search(j, cg.get(new Integer(j)));
			}
			System.out.println();
			// report.NewRound(as.getBestModelKey(j, 0), i == 0);
			// double[][] stats = as.singleBest(as.getBestModelKey(j, 0));
			// report.Add_test_statistic(stats[0]);
			// for (int k = 0; k < (i > 0 ? 0 : pr.getReplication()); k++) {
			// System.out.println(k);
			// isRabinowitzProc = true;
			// GD.RabinowitzApproach();
			// GD.CreateTableII(isRabinowitzProc);
			// GD.CreateWorkingTable(isRabinowitzProc);
			// DataFile mdrData1 = new DataFile(GD.getMarkerName(),
			// GD.getWorkingGenoTable(),
			// GD.getWorkingStatusTable(), GD.getTraitName(),
			// GD.getWorkingScoreTable(),
			// gmdrPr.getScoreIndex());
			// CombinationGenerator cg1 = new CombinationGenerator(j, j,
			// mdrData1.getMarkerNum());
			// cg1.generateCombination();
			// LinearMergeSearch as1;
			// Subdivision sd1 = new Subdivision(gmdrPr.getInterval(),
			// gmdrPr.getSeed() + i * pr.getReplication() + k,
			// mdrData1);
			// sd1.RandomPartition();
			// as1 = new LinearMergeSearch(mdrData1, sd1, cg1,
			// gmdrPr.getScoreIndex().length,
			// mdrData1.getOffset(), gmdrPr.isMooreMDR());
			// as1.search(j);
			// double[][] stats1 = as1.singleBest(as1.getBestModelKey(j, 0));
			// report.Add_null_test_statistic(stats1[0]);
			// }
			// report.RoundSummary();
			// }
			// }
			// System.out.println(report);
		}
	}
}
