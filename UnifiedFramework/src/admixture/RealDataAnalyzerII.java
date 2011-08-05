package admixture;

import java.io.IOException;
import java.util.Arrays;

import admixture.parameter.Parameter;
import admixture.parameter.ParameterParser;

import mdr.MDRConstant;
import mdr.algorithm.CombinationGenerator;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import family.mdr.AbstractMergeSearch;
import family.mdr.HeteroLinearMergeSearch;
import power.SimulationPower;
import family.pedigree.design.hierarchy.ChenInterface;
import family.pedigree.design.hierarchy.SII;
import family.pedigree.design.hierarchy.Unified;
import family.pedigree.design.hierarchy.UnifiedII;
import family.pedigree.design.hierarchy.UnifiedUnrelated;
import family.plink.PLINKParser;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class RealDataAnalyzerII {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		PLINKParser pp = new PLINKParser(p.pedigree, p.phenotype, p.map);
		long s = p.seed;
		ChenInterface chen = null;
		if (p.mode.compareTo("u") == 0) {
			if(p.unrelated_only) {
				chen = new UnifiedUnrelated(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, p.predictor, p.linkfunction);
			} else if (p.permu_fam){
				chen = new UnifiedII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, p.predictor, p.linkfunction);
			} else {
				chen = new Unified(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, p.predictor, p.linkfunction);
			}
		} else if (p.mode.compareTo("f") == 0) {
			chen = new SII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, p.predictor, p.linkfunction);
		}

		DataFile mdrData = new DataFile(chen.getMarkerName(), chen.getGenotype(), chen.getStatus(), chen.getScoreName(), chen.getScore2());
		DataFile.setScoreIndex(0);

		Subdivision sd = new Subdivision(p.cv, p.seed, mdrData.size());
		sd.RandomPartition();

		int[] includedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.includesnp);
		int[] excludedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.excludesnp);
		AbstractMergeSearch as = new HeteroLinearMergeSearch(mdrData, sd, includedMarkerIndex, excludedMarkerIndex, chen.getNumberMarker());
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
