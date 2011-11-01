package admixture.parameter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.ArrayUtils;

import family.mdr.arsenal.MDRConstant;

import test.Test;
import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Parameter {

	private final String delim = "\\s+";
	private final String incommand_separator = ",";
	private final String cmd_missing_allele = "missingallele";
	public static String missing_allele = "0";

	private final String cmd_missing_phenotype = "missingphenotype";
	public static String missing_phenotype = "-9";

	// private final String cmd_missing_genotype = "missinggenotype";
	// public static String missing_genotype = "0";

	private final String cmd_status_shift = "ss";
	public static int status_shift = 0;

	// private final String cmd_model = "model"; // cc for case control,
	// // u1 for the unified method in which the founders are exchangeable to
	// each
	// // other,
	// // u2 for the unified method in which the founders are exchangeable but
	// // within family;
	// // fam1 for ajhg2008
	// // fam2 for sii
	// public static String model = "cc";

	private final String cmd_cc = "cc";
	public static boolean ccFlag = true;

	private final String cmd_pi = "pi";
	public static boolean piFlag = false;

	private final String cmd_pii = "pii";
	public static boolean piiFlag = false;

	private final String cmd_ui = "ui";
	public static boolean uiFlag = false;

	private final String cmd_uii = "uii";
	public static boolean uiiFlag = false;

	// file set start
	private final String cmd_file = "file";
	public static boolean fileFlag = false;
	private final String cmd_ped = "ped";
	public static String ped = null;
	private final String cmd_map = "map";
	public static String map = null;

	// file set end

	// bfile set start
	private final String cmd_bfile = "bfile";
	public static boolean bfileFlag = false;
	private final String cmd_bed = "bed";
	public static String bed = null;
	private final String cmd_bim = "bim";
	public static String bim = null;
	private final String cmd_fam = "fam";
	public static String fam = null;
	// bfile set end

	// tfile set start
	// private final String cmd_tfile = "tfile";
	// public static boolean tfileFlag = false;
	// private final String cmd_tped = "tped";
	// public static String tped = null;
	// private final String cmd_tfam = "tfam";
	// public static String tfam = null;
	// tfile set end

	// phenotype set start
	private final String cmd_pheno = "pheno";
	public static String pheno = null;
	private final String cmd_res_number = "response";
	public static int response = -1;
	private final String cmd_res_name = "responsename";
	public static String response_name = null;
	private final String cmd_covar = "covar";
	public static int[] predictor = null;
	private final String cmd_covar_name = "covarname";
	public static String[] predictor_name = null;

	private final String cmd_reg = "reg";
	public int linkfunction = 0;
	// phenotype set end

	// individual selection start
	private final String cmd_ex_fam = "exfam";
	public static String[] ex_family = null;
	private final String cmd_ex_fam_file = "exfamf";
	public static boolean exfamFlag = false;
	/*
	 * private final String cmd_ex_ind = "exind"; public static String[][]
	 * ex_ind = null; private final String cmd_ex_ind_file = "exindfile"; public
	 * static boolean exindFlag = false;
	 */
	private final String cmd_filter_male = "male";
	public static boolean filter_maleFlag = false;

	private final String cmd_filter_female = "female";
	public static boolean filter_femaleFlag = false;

	private final String cmd_ex_nosex = "exnosex";
	public static boolean ex_nosexFlag = false;
	// end individual filter

	// snp selection
	private final String cmd_chr = "chr";
	public static String[] in_chr = null;
	public static String[] ex_chr = null;
	public static boolean inchrFlag = false;
	public static boolean exchrFlag = false;

	private final String cmd_snpwindow = "window";
	public static String[] snpwindow = null;
	public static double[][] snp_window = null;
	public static boolean snpwindowFlag = false;

	private final String cmd_snp = "snp";
	private final String cmd_snp_f = "snpf";
	public static boolean snpFlag = false;
	public static String[] includesnp = null;
	public static String[] excludesnp = null;

	public static String[] insnpPair = null;
	public static String[] exsnpPair = null;
	public static boolean snpPairFlag = false;

	// this set only used when the option x is specified;
	public static String[][] xincludesnp = null;
	public static String[] xinsnpPair = null;

	public static boolean xsnpFlag = false;
	public static boolean xsnpPairFlag = false;
	// end it

	private final String cmd_bgsnp = "bg";
	public static boolean bgsnpFlag = false;
	public static String[] bgsnp = null;

	private final String cmd_x = "x";
	public static boolean x = false;
	// snp selection end

	// soft snp selection
	private final String cmd_maf = "maf";
	public static double maf = 0;
	public static boolean mafFlag = false;

	private final String cmd_geno = "geno";
	public static double geno = 2;
	public static boolean genoFlag = false;

	// private final String cmd_hwe = "hwe";
	// public static double hwe = -1;
	// public static boolean hweFlag = false;

	private final String cmd_header = "header";
	public static boolean header = false;

	// sampling & partitioning start
	private final String cmd_thin = "thin";
	public static double thin = 1.0;

	private final String cmd_slice = "slice";
	public static int sliceN = 1;
	public static int slice = 1;
	public static boolean sliceFlag = false;
	// sampling & partitioning end

	// mdr options start
	private final String cmd_cv = "cv";
	public static int cv = 5;
	public static boolean cvFlag = true;

	// private final String cmd_trgroup = "trgroup";
	public static double trgroup = 0.7;
	public static boolean trgroupFlag = false;

	// private final String cmd_ttfile = "ttfile";
	public static boolean ttfileFlag = false;
	public static String[][] ttArray = null;

	// private final String cmd_trsex = "trsex";
	public static boolean trsexFlag = false;
	public static int trsex = 0;
	// private final String cmd_border = "border";
	// public static String border_fid;
	// public static String border_iid;
	// public static boolean borderFlag = false;
	// public static boolean reverseborderFlag = false;

	private final String cmd_order = "order";
	public static int order = 1;

	private final String cmd_seed = "seed";
	public static int seed = 2011;

	private final String cmd_tie = "tie";
	public static int tie = 1;
	// mdr option end

	private final String cmd_perm = "perm";
	public static int perm = 100;
	public static boolean permFlag = false;

	// private final String cmd_perm_scheme = "ps";
	public static boolean permu_scheme = true;

	private final String cmd_verbose = "verbose";
	public static boolean verboseFlag = false;

	// private final String cmd_simu = "simu";
	// public int simu = 1;

	private final String cmd_Vc = "vc";
	public static double vc = 0.05;
	public static boolean vcFlag = false;

	private final String cmd_ep = "p";
	public static double ep = 0.05;
	public static boolean epFlag = false;

	private final String cmd_training = "train";
	public static double threshold_training = 0.0;
	public static boolean trainingFlag = false;

	private final String cmd_testing = "test";
	public static double threshold_testing = 0.0;
	public static boolean testingFlag = false;

	private final String cmd_out = "out";
	public static String out = "gmdr";

	private final String cmd_help = "help";
	public boolean help = false;

	private final String cmd_testdrive = "testdrive";
	public static int testUnit = 1000;
	public static boolean testdrive = false;

	public static boolean clusterFlag = false;
	private final String cmd_node = "node";
	public static int node = 5;
	public static boolean nodeFlag = false;

	private final String cmd_email = "email";
	public static String email = "";
	public static boolean emailFlag = false;

	private final String cmd_memory = "memory";
	public static String memory = "1G";
	public static boolean memoryFlag = false;

	private final String cmd_walltime = "walltime";
	public static int walltime = 10;
	public static boolean walltimeFlag = false;

	private final String cmd_submit = "hpc";
	public static boolean submit = false;

	private final String cmd_script = "script";
	public static String script_f = "";
	private final String cmd_version = "version";
	public static String version = "\n"
			+ "************************************************************\n"
			+ "    GMDR 1.0 released 13/11/2011\n"
			+ "    Developed by Guo-Bo Chen, guobo.chen@uq.edu.au\n"
			+ "    Queensland Brain Institute, University of Queensland\n"
			+ "    St Lucia, Queensland 4067, Australia\n"
			+ "************************************************************\n";
	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();

	public Parameter() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	@SuppressWarnings("static-access")
	public void commandInitial() {
		// ops.addOption(OptionBuilder.withDescription("u (default) for the unified framework and f for using sibs only.").hasArg().create(cmd_model));
		ops.addOption(OptionBuilder.withDescription(
				"method for case-control sample.").create(cmd_cc));
		ops.addOption(OptionBuilder.withDescription(
				"method I for pedigree-based sample.").create(cmd_pi));
		ops.addOption(OptionBuilder.withDescription(
				"method II for pedigree-based sample.").create(cmd_pii));
		ops.addOption(OptionBuilder.withDescription(
				"method I for unrelated ans family samples.").create(cmd_ui));
		ops.addOption(OptionBuilder.withDescription(
				"method II for unrelated ans family samples.").create(cmd_uii));

		ops.addOption(OptionBuilder
				.withDescription("specify the .ped and .map files").hasArg()
				.create(cmd_file));
		ops.addOption(OptionBuilder.withDescription("specify the .ped file.")
				.hasArg().create(cmd_ped));
		ops.addOption(OptionBuilder.withDescription("specify the .map file.")
				.hasArg().create(cmd_map));

		ops.addOption(OptionBuilder
				.withDescription("specify the .bed, .bim and .fam files.")
				.hasArg().create(cmd_bfile));
		ops.addOption(OptionBuilder.withDescription("specify the .bed file.")
				.hasArg().create(cmd_bed));
		ops.addOption(OptionBuilder.withDescription("specify the .bim file.")
				.hasArg().create(cmd_bim));
		ops.addOption(OptionBuilder.withDescription("specify the .fam file.")
				.hasArg().create(cmd_fam));

		// ops.addOption(OptionBuilder.withDescription("specify the .tped and .tfam files").hasArg().create(cmd_tfile));
		// ops.addOption(OptionBuilder.withDescription("specify the .tped file.").hasArg().create(cmd_tped));
		// ops.addOption(OptionBuilder.withDescription("specify the .tfam file.").hasArg().create(cmd_tfam));

		ops.addOption(OptionBuilder
				.withDescription("specify the phenotype file.").hasArg()
				.create(cmd_pheno));
		ops.addOption(OptionBuilder
				.withDescription("specify 1 or more covariates by number.")
				.hasArgs().create(cmd_covar));
		ops.addOption(OptionBuilder
				.withDescription("specify 1 or more covariates by name.")
				.hasArgs().create(cmd_covar_name));

		ops.addOption(OptionBuilder
				.withDescription("specify the window size for a snp").hasArgs()
				.create(cmd_snpwindow));

		ops.addOption(OptionBuilder
				.withDescription("include snps in detecting interaction")
				.hasArgs().create(cmd_snp));
		ops.addOption(OptionBuilder
				.withDescription("specify the background snp").hasArgs()
				.create(cmd_bgsnp));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the file containing included snps when detecting interaction")
				.hasArgs().create(cmd_snp_f));
		ops.addOption(OptionBuilder.withDescription("select chromosomes")
				.hasArgs().create(cmd_chr));
		ops.addOption(OptionBuilder.withDescription("specify interacting snps")
				.create(cmd_x));

		ops.addOption(OptionBuilder
				.withDescription("specify excluded families").hasArgs()
				.create(cmd_ex_fam));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the file containing excluded family ids")
				.hasArg().create(cmd_ex_fam_file));
		// ops.addOption(OptionBuilder.withDescription("specify excluded individuals").hasArgs().create(cmd_ex_ind));
		// ops.addOption(OptionBuilder.withDescription("specify the file containing excluded individual ids").hasArg().create(cmd_ex_ind_file));
		ops.addOption(OptionBuilder.withDescription("include males only")
				.create(cmd_filter_male));
		ops.addOption(OptionBuilder.withDescription("include females only")
				.create(cmd_filter_female));
		ops.addOption(OptionBuilder.withDescription("include unknown sex ")
				.create(cmd_ex_nosex));

		ops.addOption(OptionBuilder
				.withDescription("specify response by number.").hasArg()
				.create(cmd_res_number));
		ops.addOption(OptionBuilder
				.withDescription("specify response by name.").hasArg()
				.create(cmd_res_name));

		ops.addOption(OptionBuilder
				.withDescription(
						"method for adjustment of the phenotype, 0 (default) for linear regression, 1 for logistic regression.")
				.hasArg().create(cmd_reg));
		ops.addOption(OptionBuilder
				.withDescription("fold of cross-validation, and default is 5.")
				.hasArg().create(cmd_cv));
		// ops.addOption(OptionBuilder.withDescription("specify the proportion of the training set.").hasArg().create(cmd_trgroup));
		// ops.addOption(OptionBuilder.withDescription("specify the file containing the training set.").hasArg().create(cmd_ttfile));
		// ops.addOption(OptionBuilder.withDescription("specify the gender as the training set.").hasArg().create(cmd_trsex));
		// ops.addOption(OptionBuilder.withDescription("specify start of the training set.").hasArgs().create(cmd_border));
		ops.addOption(OptionBuilder
				.withDescription("specify the order of interaction.").hasArg()
				.create(cmd_order));
		ops.addOption(OptionBuilder
				.withDescription("specify the random sample fraction.")
				.hasArg().create(cmd_thin));
		ops.addOption(OptionBuilder
				.withDescription("specify partition of the searching space.")
				.hasArg().create(cmd_slice));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the minor allele frequency for inclusion.")
				.hasArg().create(cmd_maf));
		ops.addOption(OptionBuilder
				.withDescription("specify missing genotype rate for inclusion.")
				.hasArg().create(cmd_geno));
		// ops.addOption(OptionBuilder.withDescription("specify the p value of departure from Hardy-Weinberg Equilibrium for inclusion").hasArg()
		// .create(cmd_hwe));
		ops.addOption(OptionBuilder.withDescription("seed for the algorithms")
				.hasArg().create(cmd_seed));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the classification for a tie genotype")
				.hasArg().create(cmd_tie));

		ops.addOption(OptionBuilder
				.withDescription("replication for permutation.").hasArg()
				.create(cmd_perm));
		ops.addOption(OptionBuilder
				.withDescription("output threshold for empirical p value.")
				.hasArg().create(cmd_ep));
		// ops.addOption(OptionBuilder.withDescription("only sibs are exchangeable when this option is turned on").create(cmd_perm_scheme));

		// ops.addOption(OptionBuilder.withDescription("use unrelated indivuduals only, if '--md' is specified.").create(cmd_unrelated_only));
		// ops.addOption(OptionBuilder.withDescription("replications for simulation, and this parameter is for simulation only").hasArg().create(
		// cmd_simu));
		ops.addOption(OptionBuilder
				.withDescription("missing phenotype, default -9").hasArg()
				.create(cmd_missing_phenotype));
		// ops.addOption(OptionBuilder.withDescription("missing genotype, default 00").hasArg().create(cmd_missing_genotype));
		ops.addOption(OptionBuilder
				.withDescription("missing allele, default 0").hasArg()
				.create(cmd_missing_allele));
		ops.addOption(OptionBuilder
				.withDescription(
						"use this option if status was coded as 1 (unaffected)/2 (affected).")
				.create(cmd_status_shift));
		ops.addOption(OptionBuilder
				.withDescription("threshold of variance explained for output")
				.hasArg().create(cmd_Vc));
		ops.addOption(OptionBuilder
				.withDescription("threshold of training accuracy for output")
				.hasArg().create(cmd_training));
		ops.addOption(OptionBuilder
				.withDescription("threshold of testing accuracy for output")
				.hasArg().create(cmd_testing));
		ops.addOption(OptionBuilder
				.withDescription("specify the root for output files.").hasArg()
				.create(cmd_out));
		ops.addOption(OptionBuilder.withDescription(
				"print the result in verbose form.").create(cmd_verbose));
		ops.addOption(OptionBuilder.withDescription(
				"give an evaluation for computation time")
				.create(cmd_testdrive));
		ops.addOption(OptionBuilder.withDescription("help manual.").create(
				cmd_help));
		ops.addOption(OptionBuilder.withDescription("version.").create(
				cmd_version));

		ops.addOption(OptionBuilder
				.withDescription(
						"specify the number of nodes to use on a cluster.")
				.hasArg().create(cmd_node));
		ops.addOption(OptionBuilder
				.withDescription("specify email to get informed.").hasArg()
				.create(cmd_email));
		ops.addOption(OptionBuilder.withDescription("specify memory size.")
				.hasArg().create(cmd_memory));
		ops.addOption(OptionBuilder
				.withDescription("specify wall time for each job.").hasArg()
				.create(cmd_walltime));
		ops.addOption(OptionBuilder
				.withDescription("submit jobs to a cluster.")
				.create(cmd_submit));
	}

	public void commandListenor(String[] args) {
		CommandLine cl = null;
		try {
			cl = parser.parse(ops, args);
		} catch (ParseException E) {
			System.err.println(E.getMessage());
			System.exit(0);
		}
		if (cl.hasOption(cmd_help)) {
			help = true;
		}
		if (cl.hasOption(cmd_x)) {
			x = true;
		}
		if (cl.hasOption(cmd_cc)) {
			ccFlag = true;

			piFlag = false;
			piiFlag = false;
			uiFlag = false;
			uiiFlag = false;
		}
		if (cl.hasOption(cmd_pi)) {
			piFlag = true;

			ccFlag = false;
			piiFlag = false;
			uiFlag = false;
			uiiFlag = false;
		}
		if (cl.hasOption(cmd_pii)) {
			piiFlag = true;

			ccFlag = false;
			piFlag = false;
			uiFlag = false;
			uiiFlag = false;
		}
		if (cl.hasOption(cmd_ui)) {
			uiFlag = true;

			ccFlag = false;
			piFlag = false;
			piiFlag = false;
			uiiFlag = false;
		}
		if (cl.hasOption(cmd_uii)) {
			uiiFlag = true;

			ccFlag = false;
			piFlag = false;
			piiFlag = false;
			uiFlag = false;
		}

		// file
		if (cl.hasOption(cmd_file)) {
			StringBuffer sb1 = new StringBuffer();
			StringBuffer sb2 = new StringBuffer();
			sb1.append(cl.getOptionValue(cmd_file));
			sb1.append(".ped");

			sb2.append(cl.getOptionValue(cmd_file));
			sb2.append(".map");

			ped = sb1.toString();
			map = sb2.toString();
		}
		if (cl.hasOption(cmd_ped)) {
			ped = cl.getOptionValue(cmd_ped);
		}
		if (cl.hasOption(cmd_map)) {
			map = cl.getOptionValue(cmd_map);
		}
		if (ped != null && map != null) {
			File fped = new File(ped);
			if (!fped.exists()) {
				throw new IllegalArgumentException("could not open " + ped);
			}
			File fmap = new File(map);
			if (!fmap.exists()) {
				throw new IllegalArgumentException("could not open " + map);
			}
			fileFlag = true;
		}

		// tfile
		// if (cl.hasOption(cmd_tfile)) {
		// StringBuffer sb1 = new StringBuffer();
		// StringBuffer sb2 = new StringBuffer();
		// sb1.append(cl.getOptionValue(cmd_tfile));
		// sb1.append(".tped");
		//
		// sb2.append(cl.getOptionValue(cmd_tfile));
		// sb2.append(".tfam");
		//
		// tped = sb1.toString();
		// tfam = sb2.toString();
		// }
		// if (cl.hasOption(cmd_tped)) {
		// tped = cl.getOptionValue(cmd_tped);
		// }
		// if (cl.hasOption(cmd_tfam)) {
		// tfam = cl.getOptionValue(cmd_tfam);
		// }
		// if (tped != null && tfam != null) {
		// File fped = new File(tped);
		// if (!fped.exists()) {
		// throw new IllegalArgumentException("could not open " + tped);
		// }
		// File ffam = new File(tfam);
		// if (!ffam.exists()) {
		// throw new IllegalArgumentException("could not open " + tfam);
		// }
		// tfileFlag = true;
		// }

		// bfile
		if (cl.hasOption(cmd_bfile)) {
			StringBuffer sb1 = new StringBuffer();
			StringBuffer sb2 = new StringBuffer();
			StringBuffer sb3 = new StringBuffer();
			sb1.append(cl.getOptionValue(cmd_bfile));
			sb1.append(".bed");

			sb2.append(cl.getOptionValue(cmd_bfile));
			sb2.append(".bim");

			sb3.append(cl.getOptionValue(cmd_bfile));
			sb3.append(".fam");

			bed = sb1.toString();
			bim = sb2.toString();
			fam = sb3.toString();
		}
		if (cl.hasOption(cmd_bed)) {
			bed = cl.getOptionValue(cmd_bed);
		}
		if (cl.hasOption(cmd_bim)) {
			bim = cl.getOptionValue(cmd_bim);
		}
		if (cl.hasOption(cmd_fam)) {
			fam = cl.getOptionValue(cmd_fam);
		}
		if (bed != null && bim != null && fam != null) {
			File fbed = new File(bed);
			if (!fbed.exists()) {
				throw new IllegalArgumentException("could not open " + bed);
			}
			File fbim = new File(bim);
			if (!fbim.exists()) {
				throw new IllegalArgumentException("could not open " + bim);
			}
			File ffam = new File(fam);
			if (!ffam.exists()) {
				throw new IllegalArgumentException("could not open " + fam);
			}
			bfileFlag = true;
		}

		if (cl.hasOption(cmd_pheno)) {
			pheno = cl.getOptionValue(cmd_pheno);
			File fpheno = new File(pheno);
			if (!fpheno.exists()) {
				throw new IllegalArgumentException("could not open " + fpheno);
			}
		}

		if (cl.hasOption(cmd_covar)) {
			String[] p = cl.getOptionValues(cmd_covar);
			HashSet<Integer> idx = NewIt.newHashSet();
			for (int i = 0, len = p.length; i < len; i++) {
				if (p[i].contains("-")) {
					String[] pp = p[i].split("-");
					if (pp.length != 2) {
						throw new IllegalArgumentException(
								"bad parameter for option --response ");
					}
					for (int j = Integer.parseInt(pp[0]); j <= Integer
							.parseInt(pp[1]); j++) {
						idx.add(new Integer(j));
					}
				} else {
					idx.add(new Integer(Integer.parseInt(p[i])));
				}
			}
			predictor = new int[idx.size()];
			int c = 0;
			for (Iterator<Integer> e = idx.iterator(); e.hasNext();) {
				predictor[c] = e.next().intValue() - 1;
				if (predictor[c] < 0) {
					throw new IllegalArgumentException(
							"bad parameter for option --response ");
				}
				c++;
			}
		}

		if (cl.hasOption(cmd_covar_name)) {
			HashSet<String> cn = NewIt.newHashSet();
			String[] p = cl.getOptionValues(cmd_covar_name);
			for (int i = 0; i < p.length; i++) {
				if (p[i].contains("-")) {
					String[] pp = predictor_name[i].split("-");
					if (pp.length != 2) {
						throw new IllegalArgumentException(
								"bad parameter for option --covarname ");
					}
					for (int j = 0; j < pp.length; j++) {
						cn.add(pp[j]);
					}
				} else {
					cn.add(p[i]);
				}
			}
			predictor_name = (String[]) cn.toArray(new String[0]);
		}

		if (cl.hasOption(cmd_snp)) {
			if (!x) {
				String[] snp = cl.getOptionValues(cmd_snp);
				ArrayList<String> insnp = NewIt.newArrayList();
				ArrayList<String> exsnp = NewIt.newArrayList();
				ArrayList<String> insnppair = NewIt.newArrayList();
				ArrayList<String> exsnppair = NewIt.newArrayList();
				for (int i = 0; i < snp.length; i++) {
					if (snp[i].startsWith("-")) {
						String S = snp[i].substring(1, snp[i].length());
						if (S.contains("-")) {
							String[] s = S.split("-");
							if (s.length != 2) {
								throw new IllegalArgumentException(
										"bad parameter " + snp[i]);
							}
							exsnppair.add(s[0]);
							exsnppair.add(s[1]);
						} else {
							exsnp.add(S);
						}
					} else {
						if (snp[i].contains("-")) {
							String[] s = snp[i].split("-");
							if (s.length != 2) {
								throw new IllegalArgumentException(
										"bad parameter " + snp[i]);
							}
							insnppair.add(s[0]);
							insnppair.add(s[1]);
						} else {
							insnp.add(snp[i]);
						}
					}
				}
				if (insnp.size() > 0) {
					includesnp = (String[]) insnp.toArray(new String[0]);
					snpFlag = true;
				}
				if (exsnp.size() > 0) {
					excludesnp = (String[]) exsnp.toArray(new String[0]);
					snpFlag = true;
				}
				if (insnppair.size() > 0) {
					insnpPair = (String[]) insnppair.toArray(new String[0]);
					snpPairFlag = true;
				}
				if (exsnppair.size() > 0) {
					exsnpPair = (String[]) exsnppair.toArray(new String[0]);
					snpPairFlag = true;
				}
			} else {
				String[] snp = cl.getOptionValues(cmd_snp);
				ArrayList<String> insnp = NewIt.newArrayList();
				ArrayList<String> insnppair = NewIt.newArrayList();
				for (int i = 0; i < snp.length; i++) {
					if (snp[i].startsWith("-")) {
						throw new IllegalArgumentException("bad parameter "
								+ snp[i]);
					} else {
						if (snp[i].contains("-")) {
							String[] s = snp[i].split("-");
							if (s.length != 2) {
								throw new IllegalArgumentException(
										"bad parameter " + snp[i]);
							}
							insnppair.add(s[0]);
							insnppair.add(s[1]);
						} else {
							insnp.add(snp[i]);
						}
					}
				}
				if (insnp.size() > 0) {
					xincludesnp = new String[1][];
					xincludesnp[0] = (String[]) insnp.toArray(new String[0]);
					snpFlag = true;
				}
				if (insnppair.size() > 0) {
					xinsnpPair = (String[]) insnppair.toArray(new String[0]);
					snpPairFlag = true;
				}

			}
		}

		if (cl.hasOption(cmd_bgsnp)) {
			String[] bg = cl.getOptionValues(cmd_bgsnp);
			HashSet<String> bgSet = NewIt.newHashSet();
			for (int i = 0; i < bg.length; i++) {
				bgSet.add(bg[i]);
			}
			if (bgSet.size() != bg.length) {
				throw new IllegalArgumentException("bad parameter for bgsnp");
			}
			bgsnp = cl.getOptionValues(cmd_bgsnp);
			bgsnpFlag = true;
		}

		if (cl.hasOption(cmd_snp_f)) {

			if (!x) {
				String[] snps_file = cl.getOptionValues(cmd_snp_f);
				ArrayList<String> includesnpList = NewIt.newArrayList();
				ArrayList<String> excludesnpList = NewIt.newArrayList();
				ArrayList<String> includesnpPairList = NewIt.newArrayList();
				ArrayList<String> excludesnpPairList = NewIt.newArrayList();
				for (int h = 0; h < snps_file.length; h++) {
					File f = new File(snps_file[h]);
					if (!f.exists()) {
						throw new IllegalArgumentException("could not find "
								+ snps_file[h]);
					}
					BufferedReader reader = null;
					try {
						reader = new BufferedReader(new FileReader(f));
					} catch (IOException E) {
						throw new IllegalArgumentException(
								"could not open snps file " + snps_file);
					}
					ArrayList<String> snp = NewIt.newArrayList();
					String line = null;
					try {
						while ((line = reader.readLine()) != null) {
							String[] s = line.split(delim);
							for (int i = 0; i < s.length; i++) {
								snp.add(s[i]);
							}
						}
						reader.close();
					} catch (IOException E) {
						throw new IllegalArgumentException("bad lines in "
								+ snps_file);
					}
					if (snp.size() > 0) {
						ArrayList<String> insnp = NewIt.newArrayList();
						ArrayList<String> exsnp = NewIt.newArrayList();
						ArrayList<String> insnppair = NewIt.newArrayList();
						ArrayList<String> exsnppair = NewIt.newArrayList();
						for (int i = 0; i < snp.size(); i++) {
							String subSNP = snp.get(i);
							if (subSNP.startsWith("-")) {
								String S = subSNP.substring(1, subSNP.length());
								if (S.contains("-")) {
									String[] s = S.split("-");
									if (s.length != 2) {
										throw new IllegalArgumentException(
												"bad parameter for --snpf in cmd_snp_f line "
														+ (i + 1));
									}
									exsnppair.add(s[0]);
									exsnppair.add(s[1]);
								} else {
									exsnp.add(S);
								}
							} else {
								if (subSNP.contains("-")) {
									String[] s = subSNP.split("-");
									if (s.length != 2) {
										throw new IllegalArgumentException(
												"bad parameter for --snpf in cmd_snp_f line "
														+ (i + 1));
									}
									insnppair.add(s[0]);
									insnppair.add(s[1]);
								} else {
									insnp.add(subSNP);
								}
							}
						}
						if (insnp.size() > 0) {
							includesnpList.addAll(insnp);
							snpFlag = true;
						}
						if (exsnp.size() > 0) {
							excludesnpList.addAll(exsnp);
							snpFlag = true;
						}
						if (insnppair.size() > 0) {
							includesnpPairList.addAll(insnppair);
							snpPairFlag = true;
						}
						if (exsnppair.size() > 0) {
							excludesnpPairList.addAll(exsnppair);
							snpPairFlag = true;
						}
					}
					if(includesnpList.size() >0) {
						includesnp = (String[]) includesnpList.toArray(new String[0]);
					}
					if(excludesnpList.size() > 0) {
						excludesnp = (String[]) excludesnpList.toArray(new String[0]);
					}
					if(includesnpPairList.size() > 0) {
						insnpPair = (String[]) includesnpList.toArray(new String[0]);
					}
					if(excludesnpPairList.size() > 0) {
						exsnpPair = (String[]) excludesnpList.toArray(new String[0]);
					}
				}
			} else {
				String[] snps_file = cl.getOptionValues(cmd_snp_f);
				xincludesnp = new String[snps_file.length][];
				ArrayList<String> xsnppairList = NewIt.newArrayList();
				for (int h = 0; h < snps_file.length; h++) {
					File f = new File(snps_file[h]);
					if (!f.exists()) {
						throw new IllegalArgumentException("could not find "
								+ snps_file[h]);
					}
					BufferedReader reader = null;
					try {
						reader = new BufferedReader(new FileReader(f));
					} catch (IOException E) {
						throw new IllegalArgumentException(
								"could not open snps file " + snps_file);
					}
					ArrayList<String> snp = NewIt.newArrayList();
					String line = null;
					try {
						while ((line = reader.readLine()) != null) {
							String[] s = line.split(delim);
							for (int i = 0; i < s.length; i++) {
								snp.add(s[i]);
							}
						}
						reader.close();
					} catch (IOException E) {
						throw new IllegalArgumentException("bad lines in "
								+ snps_file);
					}
					if (snp.size() > 0) {
						ArrayList<String> insnp = NewIt.newArrayList();
						ArrayList<String> insnppair = NewIt.newArrayList();

						for (int i = 0; i < snp.size(); i++) {
							String subSNP = snp.get(i);
							if (subSNP.startsWith("-")) {

							} else {
								if (subSNP.contains("-")) {
									String[] s = subSNP.split("-");
									if (s.length != 2) {
										throw new IllegalArgumentException(
												"bad parameter for --snpf in cmd_snp_f line "
														+ (i + 1));
									}
									insnppair.add(s[0]);
									insnppair.add(s[1]);
								} else {
									insnp.add(subSNP);
								}
							}
						}
						if (insnp.size() > 0) {
							xincludesnp[h] = (String[]) insnp.toArray(new String[0]);
							snpFlag = true;
						}
						if (insnppair.size() > 0) {
							xsnppairList.addAll(insnppair);
							snpPairFlag = true;
						}
					}
					if(xsnppairList.size() > 0) {
						insnpPair = (String[]) xsnppairList.toArray(new String[0]);
					}
				}
			}
		}

		if (cl.hasOption(cmd_filter_male)) {
			filter_maleFlag = true;
		}
		if (cl.hasOption(cmd_filter_female)) {
			filter_femaleFlag = true;
		}
		if (cl.hasOption(cmd_ex_nosex)) {
			ex_nosexFlag = true;
		}
		if (cl.hasOption(cmd_ex_fam)) {
			ex_family = cl.getOptionValues(cmd_ex_fam);
			HashSet<String> famSet = NewIt.newHashSet();
			for (int i = 0; i < ex_family.length; i++) {
				famSet.add(ex_family[i]);
			}
			ex_family = (String[]) famSet.toArray(new String[0]);
			exfamFlag = true;
		}
		if (cl.hasOption(cmd_ex_fam_file)) {
			String file = cl.getOptionValue(cmd_ex_fam_file);
			File f = new File(file);
			if (!f.exists()) {
				throw new IllegalArgumentException("coudl not open " + file);
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(new File(file)));
			} catch (IOException E) {
				throw new IllegalArgumentException("failed in reading " + file);
			}
			String line;
			HashSet<String> famSet = NewIt.newHashSet();
			try {
				while ((line = reader.readLine()) != null) {
					String[] l = line.split(MDRConstant.delim);
					famSet.add(l[0]);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (famSet.size() > 0) {
				ex_family = (String[]) famSet.toArray(new String[0]);
				exfamFlag = true;
			} else {
				throw new IllegalArgumentException("bad lines in " + file);
			}
		}
		/*
		 * if (cl.hasOption(cmd_ex_ind)) { String[] ind =
		 * cl.getOptionValues(cmd_ex_ind); ex_ind = new String[2][ind.length];
		 * for (int i = 0; i < ind.length / 2; i++) { String[] Ind =
		 * ind[i].split(incommand_separator); ex_ind[0][i] = Ind[0];
		 * ex_ind[1][i] = Ind[1]; } exindFlag = true; } if
		 * (cl.hasOption(cmd_ex_ind_file)) { String file =
		 * cl.getOptionValue(cmd_ex_ind_file); File find = new File(file); if
		 * (!find.exists()) { throw new
		 * IllegalArgumentException("could not open file: " + file); }
		 * ArrayList<String> fid = NewIt.newArrayList(); ArrayList<String> iid =
		 * NewIt.newArrayList();
		 * 
		 * BufferedReader reader = null; try { reader = new BufferedReader(new
		 * FileReader(new File(file))); } catch (IOException E) { throw new
		 * IllegalArgumentException("failed in reading " + file); }
		 * 
		 * String line; try { while ((line = reader.readLine()) != null) {
		 * String[] l = line.split(incommand_separator); if (l.length < 2) {
		 * continue; } fid.add(l[0]); iid.add(l[1]); } } catch (IOException e) {
		 * e.printStackTrace(System.err); System.exit(0); } if (fid.size() > 0)
		 * { exindFlag = true; ex_ind = new String[2][fid.size()]; for (int i =
		 * 0; i < ex_ind.length; i++) { ex_ind[0][i] = fid.get(i); ex_ind[1][i]
		 * = iid.get(i); } } else { throw new
		 * IllegalArgumentException("bad lines in " + file); } }
		 */
		if (cl.hasOption(cmd_filter_male)) {
			filter_maleFlag = true;
			filter_femaleFlag = false;
		}
		if (cl.hasOption(cmd_filter_female)) {
			filter_femaleFlag = true;
			filter_maleFlag = false;
		}
		if (cl.hasOption(cmd_ex_nosex)) {
			ex_nosexFlag = true;
		}

		if (cl.hasOption(cmd_chr)) {

			String[] chr = cl.getOptionValues(cmd_chr);
			HashSet<String> chrSet = NewIt.newHashSet();
			HashSet<String> exSet = NewIt.newHashSet();
			for (int i = 0; i < chr.length; i++) {
				if (chr[i].startsWith("-")) {
					exSet.add(chr[i].substring(1, chr[i].length()));
				} else {
					chrSet.add(chr[i]);
				}
			}
			if (chr.length != chrSet.size() + exSet.size()) {
				throw new IllegalArgumentException("bad parameter for --chr");
			}
			if (chrSet.size() > 0) {
				in_chr = (String[]) chrSet.toArray(new String[0]);
				inchrFlag = true;
			}
			if (exSet.size() > 0) {
				ex_chr = (String[]) exSet.toArray(new String[0]);
				exchrFlag = true;
			}
		}

		if (cl.hasOption(cmd_snpwindow)) {
			String[] s = cl.getOptionValues(cmd_snpwindow);
			snpwindow = new String[s.length];
			snp_window = new double[s.length][2];
			for (int i = 0; i < s.length; i++) {
				String[] ss = s[i].split(incommand_separator);
				if (ss.length != 3) {
					throw new IllegalArgumentException(
							"bad parameter for --snpwindow: " + s[i]);
				}
				snpwindow[i] = ss[0];
				snp_window[i][0] = Double.parseDouble(ss[1]) * -1000;
				if (Double.parseDouble(ss[2]) > 0) {
					snp_window[i][1] = Double.parseDouble(ss[2]) * 1000;
				} else {
					snp_window[i][1] = Double.MAX_VALUE;
				}
			}
			snpwindowFlag = true;
		}

		if (cl.hasOption(cmd_maf)) {
			maf = Double.parseDouble(cl.getOptionValue(cmd_maf));
			if (maf < 0) {
				throw new IllegalArgumentException("bad parameter for --maf: "
						+ maf);
			}
			mafFlag = true;
		}

		if (cl.hasOption(cmd_geno)) {
			geno = Double.parseDouble(cl.getOptionValue(cmd_geno));
			if (geno < 0) {
				throw new IllegalArgumentException("bad parameter for --geno: "
						+ geno);
			}
			genoFlag = true;
		}
		//
		// if (cl.hasOption(cmd_hwe)) {
		// hwe = Double.parseDouble(cl.getOptionValue(cmd_hwe));
		// if (hwe < 0) {
		// throw new IllegalArgumentException("bad parameter for --hwe: " +
		// hwe);
		// }
		// hweFlag = true;
		// }

		if (cl.hasOption(cmd_header)) {
			header = true;
		}
		// if (cl.hasOption(cmd_topN)) {
		// topN = Integer.parseInt(cl.getOptionValue(cmd_topN));
		// }
		if (cl.hasOption(cmd_res_number)) {
			response = Integer.parseInt(cl.getOptionValue(cmd_res_number)) - 1;
		}

		if (cl.hasOption(cmd_reg)) {
			linkfunction = Integer.parseInt(cl.getOptionValue(cmd_reg));
		}
		if (cl.hasOption(cmd_cv)) {
			cv = Integer.parseInt(cl.getOptionValue(cmd_cv));
			if (cv < 2) {
				throw new IllegalArgumentException(
						"bad parameter for option --cv.");
			}
			cvFlag = true;
		}
		// if (cl.hasOption(cmd_trgroup)) {
		// trgroup = Double.parseDouble(cl.getOptionValue(cmd_trgroup));
		// trgroupFlag = true;
		// }
		// if (cl.hasOption(cmd_trsex)) {
		// trsex = Integer.parseInt(cl.getOptionValue(cmd_trsex));
		// if (trsex != 1 && trsex != 2) {
		// throw new
		// IllegalArgumentException("unknown value for option --trsex.");
		// }
		// trsexFlag = true;
		// }
		// if (cl.hasOption(cmd_ttfile)) {
		// String tf = cl.getOptionValue(cmd_ttfile);
		// File ttfile = new File(tf);
		// if (!ttfile.exists()) {
		// throw new IllegalArgumentException("could not open ttfile " + tf);
		// }
		//
		// ArrayList<String> Farray = NewIt.newArrayList();
		// ArrayList<String> Iarray = NewIt.newArrayList();
		// BufferedReader reader = null;
		// try {
		// reader = new BufferedReader(new FileReader(new File(tf)));
		// } catch (IOException E) {
		// throw new IllegalArgumentException("failed in reading " + tf);
		// }
		// String line = null;
		// try {
		// while ((line = reader.readLine()) != null) {
		// String[] l = line.split(delim);
		// Farray.add(l[0]);
		// Iarray.add(l[1]);
		// }
		// } catch (IOException e) {
		// e.printStackTrace();
		// }
		//
		// ttArray = new String[2][Farray.size()];
		// ttArray[0] = (String[]) Farray.toArray(new String[0]);
		// ttArray[1] = (String[]) Iarray.toArray(new String[0]);
		// ttfileFlag = true;
		// }
		// if (cl.hasOption(cmd_border)) {
		// String[] h = cl.getOptionValues(cmd_border);
		// if (h.length != 2) {
		// throw new
		// IllegalArgumentException("bad parameter for option --border.");
		// }
		// boolean rflag = false;
		// if (h[0].startsWith("-")) {
		// border_fid = h[0].substring(1, h[0].length());
		// rflag = true;
		// } else {
		// border_fid = h[0];
		// }
		// if (h[1].startsWith("-")) {
		// border_iid = h[1].substring(1, h[1].length());
		// rflag = true;
		// } else {
		// border_iid = h[1];
		// }
		// borderFlag = true;
		// reverseborderFlag = rflag;
		// }
		if (cl.hasOption(cmd_seed)) {
			seed = Integer.parseInt(cl.getOptionValue(cmd_seed));
		}
		if (cl.hasOption(cmd_tie)) {
			String t = cl.getOptionValue(cmd_tie);
			if (t.compareTo("h") == 0) {
				tie = 1;
			} else if (t.compareTo("l") == 0) {
				tie = 0;
			} else {
				tie = -1;
			}
		}
		// if (cl.hasOption(cmd_simu)) {
		// simu = Integer.parseInt(cl.getOptionValue(cmd_simu));
		// }
		if (cl.hasOption(cmd_perm)) {
			perm = Integer.parseInt(cl.getOptionValue(cmd_perm));
			permFlag = true;
		}
		if (cl.hasOption(cmd_ep)) {
			ep = Double.parseDouble(cl.getOptionValue(cmd_ep));
			if (ep >= 1 || ep < 0) {
				throw new IllegalArgumentException(
						"bad parameter for option --ep.");
			}
			epFlag = true;
		}
		// if (cl.hasOption(cmd_perm_scheme)) {
		// permu_scheme = true;
		// }
		// if (cl.hasOption(cmd_unrelated_only)) {
		// unrelated_only = true;
		// }
		if (cl.hasOption(cmd_order)) {
			order = Integer.parseInt(cl.getOptionValue(cmd_order));
			if (order <= 0) {
				throw new IllegalArgumentException(
						"bad parameter for option --cv.");
			}
		}
		if (cl.hasOption(cmd_thin)) {
			thin = Double.parseDouble(cl.getOptionValue(cmd_thin));
			if (thin < 0) {
				throw new IllegalArgumentException(
						"bad parameter for option --thin.");
			}
		}
		if (cl.hasOption(cmd_slice)) {
			String[] s = cl.getOptionValue(cmd_slice).split("/");
			slice = Integer.parseInt(s[0]);
			sliceN = Integer.parseInt(s[1]);
			if (slice <= 0 || sliceN <= 0 || slice > sliceN) {
				throw new IllegalArgumentException(
						"bad parameter for option --slice.");
			}
			sliceFlag = true;
		}
		if (cl.hasOption(cmd_missing_phenotype)) {
			missing_phenotype = cl.getOptionValue(cmd_missing_phenotype);
		}
		// if (cl.hasOption(cmd_missing_genotype)) {
		// missing_genotype = cl.getOptionValue(cmd_missing_genotype);
		// }
		if (cl.hasOption(cmd_missing_allele)) {
			missing_allele = cl.getOptionValue(cmd_missing_allele);
		}
		if (cl.hasOption(cmd_status_shift)) {
			status_shift = -1;
		}
		if (cl.hasOption(cmd_Vc)) {
			vc = Double.parseDouble(cl.getOptionValue(cmd_Vc));
			vcFlag = true;
		}
		if (cl.hasOption(cmd_training)) {
			threshold_training = Double.parseDouble(cl
					.getOptionValue(cmd_training));
			trainingFlag = true;
		}
		if (cl.hasOption(cmd_testing)) {
			threshold_testing = Double.parseDouble(cl
					.getOptionValue(cmd_testing));
			testingFlag = true;
		}
		if (cl.hasOption(cmd_out)) {
			out = cl.getOptionValue(cmd_out);
		}
		if (cl.hasOption(cmd_verbose)) {
			verboseFlag = true;
		}
		if (cl.hasOption(cmd_version)) {
			System.out.println(version);
			System.exit(1);
		}
		if (cl.hasOption(cmd_testdrive)) {
			testdrive = true;
		}
		if (cl.hasOption(cmd_node)) {
			node = Integer.parseInt(cl.getOptionValue(cmd_node));
			nodeFlag = true;
			clusterFlag = true;
		}
		if (cl.hasOption(cmd_email)) {
			email = cl.getOptionValue(cmd_email);
			emailFlag = true;
			clusterFlag = true;
		}
		if (cl.hasOption(cmd_memory)) {
			memory = cl.getOptionValue(cmd_memory);
			memoryFlag = true;
			clusterFlag = true;
		}
		if (cl.hasOption(cmd_walltime)) {
			walltime = Integer.parseInt(cl.getOptionValue(cmd_walltime));
			walltimeFlag = true;
			clusterFlag = true;
		}
		if (cl.hasOption(cmd_submit)) {
			submit = true;
			clusterFlag = true;
		}
		if (help) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("UGMDR", ops);
			System.exit(1);
		}

	}

	public static void findCovar_Number(String[] F) {
		if (predictor_name != null) {
			ArrayList<Integer> idx = NewIt.newArrayList();
			for (int i = 0; i < predictor_name.length; i++) {
				if (predictor_name[i].contains("-")) {
					String[] pp = predictor_name[i].split("-");
					int j1 = ArrayUtils.indexOf(F, pp[0]);
					int j2 = ArrayUtils.indexOf(F, pp[1]);
					if (j1 > j2) {
						j1 = j1 ^ j2;
						j2 = j1 ^ j2;
						j1 = j1 ^ j2;
					}
					for (int j = j1; j <= j2; j++) {
						idx.add(new Integer(j));
					}
				} else {
					int j = ArrayUtils.indexOf(F, predictor_name[i]);
					if (j < 0) {
						System.err.println("unknown covariat "
								+ predictor_name[i]);
						System.exit(0);
					} else {
						idx.add(new Integer(j));
					}
				}
			}
			predictor = new int[idx.size()];
			for (int i = 0; i < predictor.length; i++) {
				predictor[i] = idx.get(i).intValue() - 1;
			}
		}
	}

	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}
}
