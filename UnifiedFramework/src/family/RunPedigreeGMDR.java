package family;

import java.io.IOException;
import java.util.Random;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.GMDRData;
import family.report.Report;
import mdr.data.DataFile;
import mdr.moore.LinearMergeSearch;
import mdr.GMDRParameter;

public class RunPedigreeGMDR {
	public static void main(String[] args) throws IOException {
		Report report = new Report();
		PedigreeParameter pr = new PedigreeParameter();
		long seed = 10;
		if (args.length > 0) {
			pr.read(args[0]);
		}
		GMDRParameter gmdrPr = new GMDRParameter();
		if (args.length > 1) {
			gmdrPr.read(args[1]);
		}
		boolean isRabinowitzProc = false;
		GMDRData.rnd = new Random(pr.seed);
		GMDRData GD = new GMDRData(pr.usingFounderGenotype, pr.usingChildrenGenotype, pr.isLouAlgorithm);
		AbstractGenoDistribution.rnd = new Random(pr.seed + 1);
		GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());
		if (pr.isLouAlgorithm) {
			GD.RabinowitzApproach();
		}
		GD.CreateTableII(isRabinowitzProc);
		if (pr.getScoreBuildMethod() >= 0) {
			GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(),
					pr.getScoreBuildMethod(), pr.getScoreBuildWithFounder(), pr.getScoreBuildWithChildren());
		} else {
			GD.fetchScore(pr.getPhenotypeIndex(), pr.getScoreBuildWithFounder());
		}
		GD.CreateWorkingTable(isRabinowitzProc);

		DataFile mdrData = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD.getWorkingStatusTable(), GD
				.getTraitName(), GD.getWorkingScoreTable(), gmdrPr.getScoreIndex());
		Subdivision sd = new Subdivision(gmdrPr.getInterval(), gmdrPr.getSeed(), mdrData);
		sd.RandomPartition();

		CombinationGenerator cg = new CombinationGenerator(gmdrPr.getInterctionFrom(), gmdrPr.getInteractionEnd(),
				mdrData.getMarkerNum());
		cg.generateCombination();
		LinearMergeSearch as;
		as = new LinearMergeSearch(mdrData, sd, cg, gmdrPr.getScoreIndex().length, mdrData.getOffset(), gmdrPr
				.isMooreMDR());
		for (int j = gmdrPr.getInterctionFrom(); j <= gmdrPr.getInteractionEnd(); j++) {
			as.search(j);
			report.NewRound(as.getBestModelKey(j, 0), true);
			double[][] stats = as.singleBest(as.getBestModelKey(j, 0));
			report.Add_test_statistic(stats[0]);
			for (int k = 0; k < pr.replication; k++) {
				isRabinowitzProc = true;
				GD.RabinowitzApproach();
				GD.CreateTableII(isRabinowitzProc);
				GD.CreateWorkingTable(isRabinowitzProc);
				DataFile mdrData1 = new DataFile(GD.getMarkerName(), GD.getWorkingGenoTable(), GD
						.getWorkingStatusTable(), GD.getTraitName(), GD.getWorkingScoreTable(), gmdrPr.getScoreIndex());
				CombinationGenerator cg1 = new CombinationGenerator(j, j, mdrData1.getMarkerNum());
				cg1.generateCombination();
				LinearMergeSearch as1;
				Subdivision sd1 = new Subdivision(gmdrPr.getInterval(), gmdrPr.getSeed() + k, mdrData1);
				sd1.RandomPartition();
				as1 = new LinearMergeSearch(mdrData1, sd1, cg1, gmdrPr.getScoreIndex().length, mdrData1.getOffset(),
						gmdrPr.isMooreMDR());
				as1.search(j);
				double[][] stats1 = as1.singleBest(as1.getBestModelKey(j, 0));
				report.Add_null_test_statistic(stats1[0]);
			}
			report.RoundSummary();
		}
		System.out.println(report);
	}
}
