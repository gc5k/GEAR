package admixture;

import java.io.IOException;
import java.util.Arrays;

import admixture.parameter.Parameter;
import admixture.parameter.ParameterParser;

import family.mdr.data.MDRConstant;
import family.mdr.AbstractMergeSearch;
import power.SimulationPower;
import family.mdr.HeteroCombinationSearchII;
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

public class UnifiedGMDR {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		for (int i = 0; i < p.simu; i++) {
			String PedFile = Integer.toString(i) + "L_ped.txt";
			String PhenoFile = Integer.toString(i) + "score.txt";
			String MapFile = p.map;
			PLINKParser pp = new PLINKParser(PedFile, PhenoFile, MapFile);
			
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


			int[] includedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.includesnp);
			int[] excludedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.excludesnp);
			AbstractMergeSearch as = new HeteroCombinationSearchII.Builder(Parameter.cv, chen.getSample(), chen.getMapFile(), includedMarkerIndex, excludedMarkerIndex).mute(false).build();
			for (int j = p.min; j <= p.max; j++) {

				double[] pv = new double[p.permutation];
				for (int k = 0; k < p.permutation; k++) {
					chen.getPermutedScore(p.permu_scheme);
					as.search(j, 1);
					pv[k] = as.getModelStats()[MDRConstant.TestingBalancedAccuIdx];
				}

				Arrays.sort(pv);
				chen.RecoverScore();
				as.search(j, 1);
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
