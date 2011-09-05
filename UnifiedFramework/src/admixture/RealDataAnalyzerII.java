package admixture;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import admixture.parameter.Parameter;
import admixture.parameter.ParameterParser;

import family.mdr.AbstractMergeSearch;
import family.mdr.HeteroCombinationSearchII;
import family.mdr.MDRConstant;
import power.SimulationPower;
import family.pedigree.design.hierarchy.AJHG2008;
import family.pedigree.design.hierarchy.ChenInterface;
import family.pedigree.design.hierarchy.SII;
import family.pedigree.design.hierarchy.Unified;
import family.pedigree.design.hierarchy.UnifiedII;
import family.pedigree.design.hierarchy.UnifiedUnrelated;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.plink.PLINKTransposeParser;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class RealDataAnalyzerII {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		PLINKParser pp = null;
		if(Parameter.file) {
			pp = new PLINKParser(Parameter.ped, Parameter.map, Parameter.pheno);
		} else if (Parameter.bfile){
			pp = new PLINKBinaryParser(Parameter.bed, Parameter.bim, Parameter.fam, Parameter.pheno);
		} else if (Parameter.tfile) {
			pp = new PLINKTransposeParser(Parameter.tped, Parameter.tfam, Parameter.map, Parameter.pheno);
		} else {
			System.err.println("did not specify files.");
			System.exit(0);
		}
		pp.Parse();

		long s = Parameter.seed;
		ChenInterface chen = null;
		if (Parameter.mode.compareTo("u") == 0) {
			if(p.unrelated_only) {
				chen = new UnifiedUnrelated(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, Parameter.predictor, p.linkfunction);
			} else if (p.permu_fam){
				chen = new UnifiedII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, Parameter.predictor, p.linkfunction);
			} else {
				chen = new Unified(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, Parameter.predictor, p.linkfunction);
			}
		} else if (Parameter.mode.compareTo("f") == 0) {
			chen = new SII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, Parameter.predictor, p.linkfunction);
		} else if (Parameter.mode.compareTo("pi") == 0) {
			chen = new AJHG2008(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, p.response, Parameter.predictor, p.linkfunction);
		}

		int[] includedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.includesnp);
		int[] excludedMarkerIndex = ParameterParser.selectedSNP(chen.getMapFile(), p.excludesnp);
		AbstractMergeSearch as = new HeteroCombinationSearchII.Builder(Parameter.cv, chen.getSample(), chen.getMapFile(), includedMarkerIndex, excludedMarkerIndex).mute(false).build();

		PrintStream PW = new PrintStream("ugmdr.txt");
		System.setOut(PW);
		for (int j = p.order; j <= p.order; j++) {
			System.err.println("order:" + j);
			RealDataAnalyzerII.PrintHeader(j);
//			double[] pv = new double[p.permutation];
//			for (int k = 0; k < p.permutation; k++) {
//				System.err.println("permu:" + k);
//				as.setMute(true);
//				chen.getPermutedScore(p.permu_scheme);
//				as.search(j, 1);
//				pv[k] = as.getModelStats()[MDRConstant.TestingBalancedAccuIdx];
//			}
//
//			Arrays.sort(pv);
//			double T = pv[(int) (pv.length * 0.95)];
			as.setMute(false);
			chen.RecoverScore();
			long t1 = System.currentTimeMillis();
			System.err.println(t1);
			as.search(j, 1);
			long t2 = System.currentTimeMillis();
			System.err.println(t2-t1);
//			System.out.println(as);

//			SimulationPower sp = new SimulationPower(as.getMDRResult(), pv);
//			sp.calculatePower();
//			System.out.println(sp);
		}
		System.out.println();
		PW.close();
	}

	public static void PrintHeader(int order) {
		System.out.print("The " + order + " order interaction");
		System.out.print(System.getProperty("line.separator"));
		System.out.print("model code, model(marker chr pos minor allele major allele): ");
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
