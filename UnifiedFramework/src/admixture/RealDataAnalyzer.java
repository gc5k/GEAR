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
import family.pedigree.design.hierarchy.UnifiedII;
import family.pedigree.design.hierarchy.UnifiedUnrelated;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class RealDataAnalyzer {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		String PedFile = p.pedigree;
		String PhenoFile = p.phenotype;

		ChenInterface chen = null;
		if (p.mode.compareTo("u") == 0) {
			if(p.unrelated_only) {
				chen = new UnifiedUnrelated(PedFile, PhenoFile);
			} else if (p.permu_fam){
				chen = new UnifiedII(PedFile, PhenoFile);
			} else {
				chen = new Unified(PedFile, PhenoFile);
			}
		} else if (p.mode.compareTo("f") == 0) {
			chen = new SII(PedFile, PhenoFile);
		}

		chen.setSeed(p.seed);
		chen.generateScore(p.response, p.predictor, p.linkfunction);

		DataFile mdrData = new DataFile(chen.getMarkerName(), chen.getGenotype(), chen.getStatus(), chen.getScoreName(), chen.getScore2());
		DataFile.setScoreIndex(0);

		Subdivision sd = new Subdivision(p.cv, p.seed, mdrData.size());
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
}
