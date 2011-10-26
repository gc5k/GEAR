package test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Calendar;
import java.util.List;

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
	
	public static StringBuffer LOG = new StringBuffer();
	public static void main(String[] args) throws IOException {

		Parameter p = new Parameter();
		p.commandListenor(args);
		printCommandLine(args);

		if (Parameter.clusterFlag) {
			if (!Parameter.emailFlag) {
				System.err.println("please specify your email\n");
				System.exit(1);
			}
			String script = generateScript(args);
			if (Parameter.submit) {
				Runtime   rt   =   Runtime.getRuntime();
				Process   pro   =   rt.exec(script);
				System.err.println(script + " was submitted to the cluster.");
			}
			System.exit(1);
		}
		
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
			if (Parameter.sliceFlag) {
				sb.append(".slice" + Parameter.slice + "." + Parameter.sliceN);
			}
			sb.append(".int");
			PrintStream PW = new PrintStream(sb.toString());
			System.setOut(PW);
			System.err.println("order: " + j);
			as.setMute(false);
			as.search(j, 1);
			PW.close();
			System.err.println("interaction result was saved to " + sb.toString());
			LOG.append("interaction result was saved to " + sb.toString());
			LOG.append("\n");
		}
		
		printLog();
	}
	
	public static void printCommandLine(String[] args) {
		StringBuilder sb = new StringBuilder(Parameter.out);
		sb.append(".log");

		Calendar calendar = Calendar.getInstance();
		LOG.append("The analysis was implemented at: " + calendar.getTime());
		LOG.append("\n");
		LOG.append("The command line in effect: ");
		LOG.append("\n");
		for (int i = 0; i < args.length; i++) {
			LOG.append(args[i] + " ");
		}
		LOG.append("\n");

	}
	
	public static void printLog() {
		StringBuilder sb = new StringBuilder(Parameter.out);
		if (Parameter.sliceFlag) {
			sb.append(".slice" + Parameter.slice + "." + Parameter.sliceN);
		}
		sb.append(".log");
		PrintStream pw = null;
		try {
			pw = new PrintStream(sb.toString());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		pw.append(LOG.toString());
		pw.close();
	}

	public static String generateScript(String[] args) {
		
		StringBuffer pl = new StringBuffer(Parameter.out);
		pl.append(".pl");
		PrintStream PL = null;
		try {
			PL = new PrintStream(pl.toString());
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		PL.append("#!/usr/bin/perl\n");
		PL.append("use strict;\n");
		PL.append("use warnings;\n");

		for (int i = 1; i <= Parameter.cluster; i++) {
			StringBuffer sb = new StringBuffer(Parameter.out);
			sb.append(".");
			sb.append("slice");
			sb.append(i);
			sb.append(".");
			sb.append(Parameter.cluster);
			sb.append(".");
			sb.append("cluster");
			PL.append("system(" +"\"qsub " + sb.toString() + "\"" + ");");
			PL.append("\n");
			PrintStream pw = null;
			try {
				pw = new PrintStream(sb.toString());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			pw.append("#!/bin/bash\n");
			pw.append("#$ -S /bin/bash\n");
			pw.append("#$ -cwd\n");
			pw.append("#$ -V\n");
			pw.append("#$ -m eas\n");
			pw.append("#$ -N " + Parameter.out + "." + "slice" + i + "." + Parameter.cluster + "\n");
			pw.append("#$ -M " + Parameter.email + "\n");
			pw.append("#$ -l h_rt=" + Parameter.walltime + ":10:00,s_rt=" + Parameter.walltime + ":00:00,vf=" + Parameter.memory + "\n\n");
			
			pw.append("java -jar ");
			pw.append("-Xmx" + Parameter.memory + " ");
			pw.append("gmdr.jar ");
			
			int len = args.length;
			int c = 0; 
			while(c<len){
				if(args[c].compareTo("--cluster") == 0 || args[c].compareTo("--email") == 0 || args[c].compareTo("--memory") == 0 || args[c].compareTo("--walltime") == 0) {
					c +=2;
					continue;
				}
				if(args[c].compareTo("--time") == 0) {
					c++;
					continue;
				}
				pw.append(args[c] + " ");
				c++;
			}
			pw.append("--slice ");
			pw.append(i + "/" + Parameter.cluster);
			pw.close();
			System.err.println(sb.toString() + " was generated.");
		}
		System.err.println(pl.toString() + " was generated.");
		PL.close();

		StringBuffer script = new StringBuffer();
		script.append("perl " + pl.toString());
		return script.toString();
	}
}
