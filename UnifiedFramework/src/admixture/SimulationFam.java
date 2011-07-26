package admixture;

import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

import power.SimulationPower;

import mdr.GMDRParameter;
import mdr.MDRConstant;
import mdr.algorithm.CombinationGenerator;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import mdr.moore.AbstractMergeSearch;
import mdr.moore.HeteroLinearMergeSearch;
import family.PedigreeParameter;
import family.pedigree.design.ChenAlgorithm;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class SimulationFam {

	public static void main(String[] args) throws IOException {
		int replication = args.length > 0 ? Integer.parseInt(args[0]) : 1;
		for (int i = 0; i < replication; i++) {
			String PedFile = Integer.toString(i) + "L_ped.txt";
			String PhenoFile = Integer.toString(i) + "score.txt";
			PedigreeParameter pr = new PedigreeParameter();
			long seed = 10;
			if (args.length > 1) {
				pr.read(args[1]);
				pr.setPedigreeFile(PedFile);
				pr.setPhenotypeFile(PhenoFile);
			} else {
				pr.setPedigreeFile(PedFile);
				pr.setPhenotypeFile(PhenoFile);
				pr.setCovariateIndex("3,4,5,6,7");
				pr.setPhenotypeIndex("1");
				pr.setScoreBuildMethod(1);
				pr.setAdjustScore(true);
				pr.setReplication(1);
				pr.setSeed(10);
			}
			GMDRParameter gmdrPr = new GMDRParameter();
			if (args.length > 2) {
				gmdrPr.read(args[2]);
			} else {
				gmdrPr.setInteractionEnd(2);
				gmdrPr.setInteractionFrom(2);
				gmdrPr.setInterval(10);
				gmdrPr.setSeed(seed);
				gmdrPr.setReplicationPermutation(100);
			}

			ChenAlgorithm.rnd = new Random(pr.getSeed() + i);
			ChenAlgorithm chen = new ChenAlgorithm();

			chen.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());
//			GD.print2MDRFormat("MDR.txt");
			if (pr.getScoreBuildMethod() >= 0) {
				chen.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.getScoreBuildMethod());
			} else {
				chen.fetchScore(pr.getPhenotypeIndex());
			}

			DataFile mdrData = new DataFile(chen.getMarkerName(), chen.getGenotype(), 
					chen.getStatus(), chen.getScoreName(), chen.getScore2());
			DataFile.setScoreIndex(0);
			Subdivision sd = new Subdivision(gmdrPr.getInterval(), gmdrPr.getSeed() + i, mdrData.size());
			sd.RandomPartition();

			AbstractMergeSearch as = new HeteroLinearMergeSearch(mdrData, sd);
			for (int j = gmdrPr.getInterctionFrom(); j <= gmdrPr.getInteractionEnd(); j++) {
				CombinationGenerator cg = new CombinationGenerator(j, mdrData.getNumMarker());
				cg.generateCombination();

				double[] p = new double[gmdrPr.getReplicationPermutation()];
				for (int k = 0; k < gmdrPr.getReplicationPermutation(); k++) {
					mdrData.setScore(chen.getPermutedScore(true));
					as.search(j, cg.getCombination());
					p[k] = as.getStats()[MDRConstant.TestingBalancedAccuIdx];
				}

				Arrays.sort(p);
				mdrData.setScore(chen.getScore());
				as.search(j, cg.getCombination());
				System.out.println(as);

				SimulationPower sp = new SimulationPower(as.getMDRResult(), p);
				sp.calculatePower();
				System.out.println(sp);
			}
		}
	}
}
