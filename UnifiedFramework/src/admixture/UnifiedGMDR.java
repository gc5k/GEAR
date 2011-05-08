package admixture;

import java.io.IOException;
import java.util.Arrays;

import admixture.parameter.Parameter;

import mdr.MDRConstant;
import mdr.algorithm.CombinationGenerator;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import mdr.moore.AbstractMergeSearch;
import mdr.moore.HeteroLinearMergeSearch;
import power.SimulationPower;
import family.pedigree.design.hierarchy.ChenInterface;
import family.pedigree.design.hierarchy.SII;
import family.pedigree.design.hierarchy.Unified;

public class UnifiedGMDR {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		for (int i = 0; i < p.simu; i++) {
			String PedFile = Integer.toString(i) + "L_ped.txt";
			String PhenoFile = Integer.toString(i) + "score.txt";

			ChenInterface chen = null;
			if (p.mode.compareTo("u") == 0) {
				chen = new Unified(PedFile, PhenoFile);
			} else if (p.mode.compareTo("f") == 0) {
				chen = new SII(PedFile, PhenoFile);
			}
			chen.setSeed(p.seed);
			chen.generateScore(p.response, p.predictor, p.linkfunction);

			DataFile mdrData = new DataFile(chen.getMarkerName(), chen.getGenotype(), chen.getStatus(), chen.getScoreName(), chen.getScore2());
			DataFile.setScoreIndex(0);

			Subdivision sd = new Subdivision(p.cv, p.seed + i, mdrData.size());
			sd.RandomPartition();

			AbstractMergeSearch as = new HeteroLinearMergeSearch(mdrData, sd);
			for (int j = p.min; j <= p.max; j++) {
				CombinationGenerator cg = new CombinationGenerator(j, mdrData.getNumMarker());
				cg.generateCombination();

				double[] pv = new double[p.permutation];
				for (int k = 0; k < p.permutation; k++) {
					mdrData.setScore(chen.getPermutedScore(p.permu_scheme));
					as.search(j, cg.getCombination());
					pv[k] = as.getStats()[MDRConstant.TestingBalancedAccuIdx];
				}

				Arrays.sort(pv);
				mdrData.setScore(chen.getScore());
				as.search(j, cg.getCombination());
				System.out.println(as);

				SimulationPower sp = new SimulationPower(as.getMDRResult(), pv);
				sp.calculatePower();
				System.out.println(sp);
			}
		}
		System.out.println("type I error at alpha=0.05: " + SimulationPower.typeI_005 + "/" + p.simu);
		System.out.println("type I error at alpha=0.01: " + SimulationPower.typeI_001 + "/" + p.simu);
		System.out.println("power at alpha = 0.05: " + SimulationPower.power_005 + "/" + p.simu);
		System.out.println("power at alpha = 0.01: " + SimulationPower.power_001 + "/" + p.simu);
	}
}
