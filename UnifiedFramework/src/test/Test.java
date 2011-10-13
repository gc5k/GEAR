package test;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import admixture.parameter.Parameter;

import family.mdr.AbstractMergeSearch;
import family.mdr.HeteroCombinationSearchII;
import family.mdr.arsenal.MDRConstant;
import family.mdr.arsenal.ModelGenerator;
import family.mdr.arsenal.ModelGeneratorII;
import family.mdr.filter.softfilter.SoftSNPFilter;
import family.pedigree.design.hierarchy.AJHG2008;
import family.pedigree.design.hierarchy.ChenInterface;
import family.pedigree.design.hierarchy.SII;
import family.pedigree.design.hierarchy.Unified;
import family.pedigree.design.hierarchy.UnifiedII;
import family.pedigree.design.hierarchy.UnifiedUnrelated;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.plink.PLINKTransposeParser;
import family.popstat.AlleleFrequency;
import family.popstat.GenotypeMatrix;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Test {
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);

		PLINKParser pp = null;
		if (Parameter.fileFlag) {
			pp = new PLINKParser(Parameter.ped, Parameter.map, Parameter.pheno);
		} else if (Parameter.bfileFlag) {
			pp = new PLINKBinaryParser(Parameter.bed, Parameter.bim, Parameter.fam, Parameter.pheno);
		} else if (Parameter.tfileFlag) {
			pp = new PLINKTransposeParser(Parameter.tped, Parameter.tfam, Parameter.map, Parameter.pheno);
		} else {
			System.err.println("did not specify files.");
			System.exit(0);
		}
		pp.Parse();

		long s = Parameter.seed;
		ChenInterface chen = null;
		if (Parameter.model.compareTo("cc") == 0) {
			chen = new UnifiedUnrelated(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor,
						p.linkfunction);
		} else if (Parameter.model.compareTo("u1") == 0) {
			chen = new Unified(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor, p.linkfunction);
		} else if (Parameter.model.compareTo("u2") == 0) {
			chen = new UnifiedII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor, p.linkfunction);	
		} else if (Parameter.model.compareTo("fam1") == 0) {
			chen = new AJHG2008(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor, p.linkfunction);
		} else if (Parameter.model.compareTo("fam2") == 0) {
			chen = new SII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor, p.linkfunction);
		}

		GenotypeMatrix gm = new GenotypeMatrix(chen);
//		gm.Test();
		AlleleFrequency AF = new AlleleFrequency(gm);
		AF.CalculateAlleleFrequency();
		System.out.println(AF);

		GenotypeMatrix GM = new GenotypeMatrix(chen);
		AlleleFrequency af = new AlleleFrequency(GM);
		af.CalculateAlleleFrequency();
		pp.setAlleleFrequency(af.getAlleleFrequency());

		SoftSNPFilter snpFilterII = new SoftSNPFilter(pp.getSNPFilter(), af);

		AbstractMergeSearch as;
		ModelGenerator mg;
		if (Parameter.x) {
			mg = new ModelGeneratorII(snpFilterII.getWSeq2(), snpFilterII.getBgSeq());
		} else {
			mg = new ModelGenerator(snpFilterII.getWSeq(), snpFilterII.getBgSeq());
		}
		as = new HeteroCombinationSearchII.Builder(Parameter.cv, chen.getSample(), chen.getMapFile()).
		ModelGenerator(mg).mute(false).build();
		
//		PrintStream PW = new PrintStream("ugmdr.txt");
//		System.setOut(PW);
//		for (int j = p.order; j <= p.order; j++) {
//			System.err.println("order:" + j);
//			Test.PrintHeader(j);
//			// double[] pv = new double[p.permutation];
//			// for (int k = 0; k < p.permutation; k++) {
//			// System.err.println("permu:" + k);
//			// as.setMute(true);
//			// chen.getPermutedScore(p.permu_scheme);
//			// as.search(j, 1);
//			// pv[k] = as.getModelStats()[MDRConstant.TestingBalancedAccuIdx];
//			// }
//			//
//			// Arrays.sort(pv);
//			// double T = pv[(int) (pv.length * 0.95)];
//			as.setMute(false);
//			chen.RecoverScore();
//			long t1 = System.currentTimeMillis();
//			System.err.println(t1);
//			as.search(j, 1);
//			long t2 = System.currentTimeMillis();
//			System.err.println(t2 - t1);
//			// System.out.println(as);
//
//			// SimulationPower sp = new SimulationPower(as.getMDRResult(), pv);
//			// sp.calculatePower();
//			// System.out.println(sp);
//		}
//		System.out.println();
//		PW.close();
	}

	public static void PrintHeader(int order) {
		System.out.print("The " + order + " order interaction");
		System.out.print(System.getProperty("line.separator"));
		System.out.print("model code, model(marker chr pos minor allele major allele): ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if (i != MDRConstant.NumStats - 1) {
				System.out.print(MDRConstant.TestStatistic[i] + ", ");
			} else {
				System.out.print(MDRConstant.TestStatistic[i]);
			}
		}
		System.out
				.print(": classfication (genotype, High-risk or Low-risk group, positive scores, positive subjects, negative score, negative subjects)");
		System.out.println();
	}
}
