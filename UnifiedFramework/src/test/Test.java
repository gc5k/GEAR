package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Calendar;

import util.NewIt;

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
import family.plink.PLINKParser; //import family.plink.PLINKTransposeParser;
import family.popstat.AlleleFrequency;
import family.popstat.GenotypeMatrix;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Test {

	public static StringBuffer LOG = new StringBuffer();
	public static boolean fileFlag = false;
	public static boolean bfileFlag = false;
	public static void main(String[] args) throws IOException {

		System.err.print(Parameter.version);
		String[] scmd = script(args);
		Parameter p = new Parameter();
		printCommandLine(scmd);
		p.commandListenor(scmd);

//		savecmd(args);
		if (Parameter.clusterFlag) {
			String script = generateScript(scmd);

			if (Parameter.submit) {
				Runtime rt = Runtime.getRuntime();
				rt.exec(script);
				System.err.println(script + " was submitted to the cluster.");
				Test.LOG.append(script + " was submitted to the cluster.\n");
			}
			Test.printLog();
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
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
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
		// System.err.println(af);
		pp.setAlleleFrequency(af.getAlleleFrequency());

		SoftSNPFilter softFilter = new SoftSNPFilter(pp.getSNPFilter(), af);
		softFilter.Filter();

		AbstractMergeSearch as;
		ModelGenerator mg;
		if (Parameter.transFlag) {
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
			// sb.append(j);
			if (Parameter.sliceFlag) {
				sb.append(".slice" + Parameter.slice + "." + Parameter.sliceN);
			}
			sb.append(".int");
			PrintStream PW = new PrintStream(sb.toString());
			System.setOut(PW);
			System.err.println("interaction order: " + j);
			as.setMute(false);
			as.search(j, 1);
			PW.close();
			System.err.println("\ninteraction result was saved to " + sb.toString());
			LOG.append("\ninteraction result was saved to " + sb.toString() + "\n");

		}

		printLog();
	}

	public static void printCommandLine(String[] args) {

		Calendar calendar = Calendar.getInstance();
		System.err.println("The analysis was implemented at: " + calendar.getTime() + "\n");
		LOG.append(Parameter.version);
		LOG.append("The analysis was implemented at: " + calendar.getTime() + "\n\n");
		System.err.println("The command line in effect: ");
		LOG.append("The command line in effect: \n");
		for (int i = 0; i < args.length; i++) {
			System.err.print(args[i] + " ");
			LOG.append(args[i] + " ");
		}
		System.err.println("\n");
		LOG.append("\n\n");

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
		Calendar calendar = Calendar.getInstance();
		System.err.println("\nThe analysis was finished at: " + calendar.getTime() + "\n");
		System.err.println("These above messages were printed into " + sb.toString() + ".\n");
		LOG.append("\nThe analysis was finished at: " + calendar.getTime() + "\n");
		LOG.append("\n");

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

		for (int i = 1; i <= Parameter.node; i++) {
			StringBuffer sb = new StringBuffer(Parameter.out);
			sb.append(".");
			sb.append("slice");
			sb.append(i);
			sb.append(".");
			sb.append(Parameter.node);
			sb.append(".");
			sb.append("cluster");
			PL.append("system(" + "\"qsub " + sb.toString() + "\"" + ");");
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
			pw.append("#$ -N " + Parameter.out + "." + "slice" + i + "." + Parameter.node + "\n");
			if (Parameter.emailFlag) {
				pw.append("#$ -M " + Parameter.email + "\n");
			}
			pw.append("#$ -l h_rt=" + Parameter.walltime + ":10:00,s_rt=" + Parameter.walltime + ":00:00,vf=" + Parameter.memory + "\n\n");

			pw.append("java -jar ");
			pw.append("-Xmx" + Parameter.memory + " ");
			pw.append("gmdr.jar ");

			int len = args.length;
			int c = 0;
			while (c < len) {
				if (args[c].compareTo("--node") == 0 || args[c].compareTo("--email") == 0 || args[c].compareTo("--memory") == 0
						|| args[c].compareTo("--walltime") == 0) {
					c += 2;
					continue;
				}
				if (args[c].compareTo("--testdrive") == 0) {
					c++;
					continue;
				}
				if (args[c].compareTo("--hpc") == 0) {
					c++;
					continue;
				}
				pw.append(args[c] + " ");
				c++;
			}
			pw.append("--slice ");
			pw.append(i + "/" + Parameter.node);
			pw.close();
			System.err.println(sb.toString() + " was generated.");
		}
		System.err.println(pl.toString() + " was generated.");
		PL.close();

		StringBuffer script = new StringBuffer();
		script.append("perl " + pl.toString());
		return script.toString();
	}

	public static String[] script(String[] args) {

		String sf = null;
		boolean scriptFlag = false;
		String[] scmd;
		int c = 0;
		for (int i = 0; i < args.length; i++) {
			if (args[i].compareTo("--script") == 0) {
				c = i;
				sf = args[i + 1];
				scriptFlag = true;
				continue;
			} else if (i == c+1) {
				continue;
			}
			if (args[i].compareTo("--bfile")==0 || args[i].compareTo("--bed") == 0 || args[i].compareTo("--bim") == 0 || args[i].compareTo("--fam") == 0) {
				bfileFlag = true;
			}
			if (args[i].compareTo("--file") == 0 || args[i].compareTo("--ped") == 0 || args[i].compareTo("--map") == 0) {
				fileFlag = true;
			}
		}

		if (bfileFlag && fileFlag) {
			System.err.println("specified both text and binary format files.");
			Test.LOG.append("specified both text and binary format files.\n");
			Test.printLog();
			System.exit(0);
		}

		if (scriptFlag) {
			File f = new File(sf);
			if (!f.exists()) {
				System.err.print("could not find " + sf + ".");
				Test.LOG.append("could not find " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(f));
			} catch (IOException E) {
				System.err.print("could not open " + sf + ".");
				Test.LOG.append("could not open " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}

			ArrayList<String> cmd = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split("\\s+");
					for (int i = 0; i < s.length; i++) {
						cmd.add(s[i]);
					}
				}
				reader.close();
			} catch (IOException E) {
				System.err.println("bad lines in " + sf + ".");
				Test.LOG.append("bad lines in " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}

			scmd = (String[]) cmd.toArray(new String[0]);
		} else {
			scmd = args;
		}
		return scmd;
	}

	public static void savecmd(String[] args) {
		StringBuffer sb = new StringBuffer();

		sb.append(Parameter.out);
		if(Parameter.submit) {
			sb.append(".hpc");
		} else if(Parameter.clusterFlag) {
			sb.append(".");
			sb.append("slice");
			sb.append(Parameter.slice);
			sb.append(".");
			sb.append(Parameter.sliceN);
			sb.append(".");
			sb.append("cluster");
		}

		sb.append(".script");
		PrintStream ps = null;
		try {
			ps = new PrintStream(sb.toString());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		for (String s:args) {
			ps.print(s + " ");
		}
		ps.close();
	}
}
