package admixture;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import admixture.parameter.Parameter;
import admixture.parameter.ParameterParser;

import mdr.MDRConstant;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import family.mdr.AbstractMergeSearch;
import family.mdr.HeteroCombinationSearchII;
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
		AbstractMergeSearch as = new HeteroCombinationSearchII.Builder(mdrData, sd, chen.getMapFile(), includedMarkerIndex, excludedMarkerIndex).mute(false).build();

		PrintStream PW = new PrintStream("ugmdr.txt");
		System.setOut(PW);
		for (int j = p.min; j <= p.max; j++) {
			RealDataAnalyzerII.PrintHeader(j);
			double[] pv = new double[p.permutation];
			for (int k = 0; k < p.permutation; k++) {
				as.setMute(true);
				mdrData.setScore(chen.getPermutedScore(p.permu_scheme));
				as.search(j, 1);
				pv[k] = as.getModelStats()[MDRConstant.TestingBalancedAccuIdx];
			}
			
			Arrays.sort(pv);
			double T = pv[(int) (pv.length * 0.95)];
			as.setMute(false);
			mdrData.setScore(chen.getScore());
			as.search(j, 1);
//			System.out.println(as);

			SimulationPower sp = new SimulationPower(as.getMDRResult(), pv);
			sp.calculatePower();
			System.out.println(sp);
		}
		System.out.println();
		PW.close();
	}

	public static void PrintHeader(int order) {
		System.out.print("The " + order + " order interaction");
		System.out.print(System.getProperty("line.separator"));
		System.out.print("model(marker chr pos): ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if( i != MDRConstant.NumStats - 1) {
				System.out.print(MDRConstant.TestStatistic[i] + ", ");
			} else {
				System.out.print(MDRConstant.TestStatistic[i]);
			}
		}
		System.out.print(": classfication (genotype, High-risk or Low-risk group, positive scores, positive subjects, negative score, negative subjects)");
		System.out.println();
	}
}
