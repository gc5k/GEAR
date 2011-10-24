package test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Calendar;
import admixture.parameter.Parameter;

import family.mdr.AbstractMergeSearch;
import family.mdr.HeteroCombinationSearchP;
import family.mdr.TTMDR;
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
		printCommandLine(args);

		PLINKParser pp = null;
		if (Parameter.fileFlag) {
			pp = new PLINKParser(Parameter.ped, Parameter.map, Parameter.pheno);
		} else if (Parameter.bfileFlag) {
			pp = new PLINKBinaryParser(Parameter.bed, Parameter.bim, Parameter.fam, Parameter.pheno);
			// } else if (Parameter.tfileFlag) {
			// pp = new PLINKTransposeParser(Parameter.tped, Parameter.tfam,
			// Parameter.map, Parameter.pheno);
		} else {
			System.err.println("did not specify files.");
			System.exit(0);
		}
		pp.Parse();

		long s = Parameter.seed;
		ChenInterface chen = null;
		if (Parameter.ccFlag) {
			chen = new UnifiedUnrelated(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor,
					p.linkfunction);
		} else if (Parameter.uiFlag) {
			chen = new Unified(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor,
					p.linkfunction);
		} else if (Parameter.uiiFlag) {
			chen = new UnifiedII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor,
					p.linkfunction);
		} else if (Parameter.piFlag) {
			chen = new AJHG2008(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor,
					p.linkfunction);
		} else if (Parameter.piiFlag) {
			chen = new SII(pp.getPedigreeData(), pp.getPhenotypeData(), pp.getMapData(), s, Parameter.response, Parameter.predictor, p.linkfunction);
		}

		GenotypeMatrix GM = new GenotypeMatrix(chen);
		AlleleFrequency af = new AlleleFrequency(GM);
		af.CalculateAlleleFrequency();
//		System.err.println(af);
		pp.setAlleleFrequency(af.getAlleleFrequency());

		SoftSNPFilter softFilter = new SoftSNPFilter(pp.getSNPFilter(), af);
		softFilter.Filter();

		AbstractMergeSearch as;
		ModelGenerator mg;
		if (Parameter.x) {
			mg = new ModelGeneratorII(softFilter.getWSeq2(), softFilter.getBgSeq());
		} else {
			mg = new ModelGenerator(softFilter.getWSeq(), softFilter.getBgSeq());
		}
		if (Parameter.trgroupFlag) {
			as = new TTMDR(2, chen.getSample(), chen.getMapFile(), mg, 1, false);
		} else {
			as = new HeteroCombinationSearchP.Builder(Parameter.cv, chen.getSample(), chen.getMapFile()).ModelGenerator(mg).mute(false).chen(chen)
					.build();
		}

		for (int j = Parameter.order; j <= Parameter.order; j++) {
			StringBuilder sb = new StringBuilder(Parameter.out);
//			sb.append(j);
			sb.append(".int");
			if (Parameter.sliceFlag) {
				sb.append(".slice" + Parameter.slice + "." + Parameter.sliceN);
			}
			PrintStream PW = new PrintStream(sb.toString());
			System.setOut(PW);
			System.err.println("order: " + j);
			as.setMute(false);
			as.search(j, 1);
			PW.close();
			System.out.println("interaction result was saved to " + sb.toString());
			System.err.println(System.currentTimeMillis());
		}
	}
	
	public static void printCommandLine(String[] args) {
		StringBuilder sb = new StringBuilder(Parameter.out);
		sb.append(".log");

		PrintStream PW = null;
		try {
			PW = new PrintStream(sb.toString());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Calendar calendar = Calendar.getInstance();
		PW.println("The analysis was implemented at: " + calendar.getTime());
		PW.println("The command line in effect:");
		for (int i = 0; i < args.length; i++) {
			PW.print(args[i] + " ");
		}
		PW.close();
	}
}
