package admixture.parameter;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.GnuParser;
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
	private final String cmd_missing_allele_long = "missing-allele";
	public static String missing_allele = "0";

	private final String cmd_missing_phenotype = "missingphenotype";
	private final String cmd_missing_phenotype_long = "missing-phenotype";
	public static String missing_phenotype = "-9";

	// private final String cmd_missing_genotype = "missinggenotype";
	// public static String missing_genotype = "0";

	private final String cmd_status_shift = "1";
	public static int status_shift = 1;
	public static boolean status_shiftFlag = false;

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
	private final String cmd_covar = "covar";
	public static String pheno = null;
	private final String cmd_pheno_number = "pheno_number";
	private final String cmd_pheno_number_long = "pheno-number";
	public static int response = -1;
	private final String cmd_pheno_name = "response_name";
	private final String cmd_pheno_name_long = "pheno-name";
	public static String response_name = null;
	private final String cmd_covar_number = "covar_number";
	private final String cmd_covar_number_long = "covar-number";
	public static int[] predictor = null;
	private final String cmd_covar_name = "covar_name";
	private final String cmd_covar_name_long = "covar-name";
	public static String[] predictor_name = null;

	private final String cmd_reg = "regression";
	public int linkfunction = 0;
	// phenotype set end

	// individual selection start
	private final String cmd_remove = "remove";
	public static String[][] ex_family = null;
	public static boolean removeFlag = false;
	private final String cmd_keep = "keep";
	public static String[][] indKeep = null;
	public static boolean keepFlag = false;
		
	
	/*
	 * private final String cmd_ex_ind = "exind"; public static String[][]
	 * ex_ind = null; private final String cmd_ex_ind_file = "exindfile"; public
	 * static boolean exindFlag = false;
	 */
	private final String cmd_keep_male = "male";
	private final String cmd_keep_male_long = "keep-male";
	public static boolean keep_maleFlag = false;

	private final String cmd_keep_female = "female";
	private final String cmd_keep_female_long = "keep-female";
	public static boolean keep_femaleFlag = false;

	private final String cmd_ex_nosex = "exnosex";
	private final String cmd_ex_nosex_long = "exclude-nosex";
	public static boolean ex_nosexFlag = false;
	// end individual filter

	// snp selection
	private final String cmd_region = "region";
	public static boolean regionFlag = false;
	public static String[] chr_reg = null;
	public static double[] begin = null;
	public static double[] end = null;

	private final String cmd_gene_window = "genewindow";
	private final String cmd_gene_window_long = "gene-window";
	public static double genewindow = 0;
	
	private final String cmd_hg18 = "hg18";
	public static boolean hg18Flag = false;
	private final String cmd_hg19 = "hg19";
	public static boolean hg19Flag = true;
	private final String cmd_hg = "hg";
	public static String hgFile = "/gene36.txt";
	
	private final String cmd_gene = "gene";
	private final String cmd_gene_list = "genelist";
	private final String cmd_gene_list_long = "gene-list";
	public static boolean geneFlag = false;
	public static boolean genelistFlag = false;
	public static String[] gene = null;
	public static String[] gene_chr = null;
	public static double[] gene_begin = null;
	public static double[] gene_end = null;

	private final String cmd_snp2genelist = "makegene2snplist";
	private final String cmd_snp2gene_list = "make-snp-list";
	private final String cmd_snp2genemlist = "makegene2snpmlist";
	private final String cmd_snp2gene_mlist = "make-snp-mlist";
	public static boolean snp2genefilesFlag = false;
	public static boolean snp2genefileFlag = false;
	
	private final String cmd_chr = "chr";
	public static String[] in_chr = null;
	public static String[] ex_chr = null;
	public static boolean inchrFlag = false;
	public static boolean exchrFlag = false;

	private final String cmd_snpwindow = "snpwindow";
	private final String cmd_snpwindow_long = "snp-window";
	public static String[] snpwindow = null;
	public static double[][] snp_window = null;
	public static boolean snpwindowFlag = false;

//	private final String cmd_snp = "snp";
	private final String cmd_extract = "extract";
	public static boolean snpFlag = false;
	private final String cmd_exclude = "exclude";
	public static String[] includesnp = null;
	public static String[] excludesnp = null;

	// this set only used when the option x is specified;
	public static String[][] xincludesnp = null;
	public static boolean xsnpFlag = false;
	// end it

	private final String cmd_bgsnp = "bg";
	public static boolean bgsnpFlag = false;
	public static String[] bgsnp = null;

	private final String cmd_trans = "trans";
	public static boolean transFlag = false;
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
	private final String cmd_testdrive_long = "test-drive";
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

	public static String script_f = "";
	private final String cmd_version = "version";
	public static String version = "\n"
			+ "******************************************************************\n"
			+ "| GMDR 1.0 released 13/11/2011                                   |\n"
			+ "| (C) 2011 Guo-Bo Chen, Xiang-Yang Lou                           |\n"
			+ "| v 0.7.3                                                        |\n"			
			+ "| GNU General Public License, v2                                 |\n"
			+ "| Department of Biostatistics, Section on Statistical Genetics   |\n"
			+ "| University of Alabama at Birmingham                            |\n"
			+ "******************************************************************\n";
	private Options ops = new Options();
	private CommandLineParser parser = new GnuParser();

	public Parameter() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	@SuppressWarnings("static-access")
	public void commandInitial() {

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

		ops.addOption(OptionBuilder
				.withDescription("specify the phenotype file.").hasArg()
				.create(cmd_covar));
		ops.addOption(OptionBuilder
				.withDescription("specify 1 or more covariates by number.")
				.hasArgs().withLongOpt(cmd_covar_number_long).create(cmd_covar_number));
		ops.addOption(OptionBuilder
				.withDescription("specify 1 or more covariates by name.")
				.hasArgs().withLongOpt(cmd_covar_name_long).create(cmd_covar_name));

		ops.addOption(OptionBuilder
				.withDescription("specify regions to select snps.").hasArgs()
				.create(cmd_region));
		ops.addOption(OptionBuilder.withDescription("specify genes.").hasArgs()
				.create(cmd_gene));
		ops.addOption(OptionBuilder.withDescription("specify a genelist.").hasArgs()
				.withLongOpt(cmd_gene_list_long).create(cmd_gene_list));
		ops.addOption(OptionBuilder.withDescription("specify gene window in kb.").hasArg()
				.withLongOpt(cmd_gene_window_long).create(cmd_gene_window));

		ops.addOption(OptionBuilder.withDescription("make snp lists with respect genes")
				.withLongOpt(cmd_snp2gene_list).create(cmd_snp2genelist));
		ops.addOption(OptionBuilder.withDescription("make snp lists with respect genes")
				.withLongOpt(cmd_snp2gene_mlist).create(cmd_snp2genemlist));

		ops.addOption(OptionBuilder.withDescription("specify human genome.").hasArg()
				.create(cmd_hg));
		ops.addOption(OptionBuilder.withDescription("use human genome build 18.")
				.create(cmd_hg18));
		ops.addOption(OptionBuilder.withDescription("use human genome build 19.")
				.create(cmd_hg19));
		

		ops.addOption(OptionBuilder
				.withDescription("specify the window size for a snp").hasArgs()
				.withLongOpt(cmd_snpwindow_long).create(cmd_snpwindow));
		ops.addOption(OptionBuilder
				.withDescription("specify the background snp").hasArgs()
				.create(cmd_bgsnp));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the file containing included snps when detecting interaction")
				.hasArgs().create(cmd_extract));
		ops.addOption(OptionBuilder
				.withDescription(
						"specify the file containing excluded snps when detecting interaction")
				.hasArg().create(cmd_exclude));
		ops.addOption(OptionBuilder.withDescription("select chromosomes")
				.hasArgs().create(cmd_chr));
		ops.addOption(OptionBuilder.withDescription("specify interacting snps")
				.create(cmd_trans));

		ops.addOption(OptionBuilder
				.withDescription("remove individuals").hasArg()
				.create(cmd_remove));

		ops.addOption(OptionBuilder
				.withDescription("keep individuals").hasArg()
				.create(cmd_keep));
		// ops.addOption(OptionBuilder.withDescription("specify excluded individuals").hasArgs().create(cmd_ex_ind));
		// ops.addOption(OptionBuilder.withDescription("specify the file containing excluded individual ids").hasArg().create(cmd_ex_ind_file));
		ops.addOption(OptionBuilder.withDescription("keep males only").withLongOpt(cmd_keep_male_long)
				.create(cmd_keep_male));
		ops.addOption(OptionBuilder.withDescription("keep females only").withLongOpt(cmd_keep_female_long)
				.create(cmd_keep_female));
		ops.addOption(OptionBuilder.withDescription("exclude unknown sex").withLongOpt(cmd_ex_nosex_long)
				.create(cmd_ex_nosex));

		ops.addOption(OptionBuilder
				.withDescription("specify phenotype by number.").hasArg()
				.withLongOpt(cmd_pheno_number_long).create(cmd_pheno_number));
		ops.addOption(OptionBuilder
				.withDescription("specify phenotype by name.").hasArg()
				.withLongOpt(cmd_pheno_name_long).create(cmd_pheno_name));

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
				.withDescription("missing phenotype, default -9").hasArg().withLongOpt(cmd_missing_phenotype_long)
				.create(cmd_missing_phenotype));
		// ops.addOption(OptionBuilder.withDescription("missing genotype, default 00").hasArg().create(cmd_missing_genotype));
		ops.addOption(OptionBuilder
				.withDescription("missing allele, default 0").hasArg().withLongOpt(cmd_missing_allele_long)
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
				"give an evaluation for computation time").withLongOpt(cmd_testdrive_long)
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
		if (cl.hasOption(cmd_trans)) {
			transFlag = true;
		}
		if (cl.hasOption(cmd_out)) {
			out = cl.getOptionValue(cmd_out);
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
				System.err.println("could not open " + ped + ".");
				Test.LOG.append("could not open " + ped + ".\n");
				Test.printLog();
				System.exit(0);
			}
			File fmap = new File(map);
			if (!fmap.exists()) {
				System.err.println("could not open " + map + ".");
				Test.LOG.append("could not open " + map + ".\n");
				Test.printLog();
				System.exit(0);
			}
			fileFlag = true;
		}

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
				System.err.println("could not open " + bed + ".");
				Test.LOG.append("could not open " + bed + ".\n");
				Test.printLog();
				System.exit(0);
			}
			File fbim = new File(bim);
			if (!fbim.exists()) {
				System.err.println("could not open " + bim +".");
				Test.LOG.append("could not open " + bim + ".\n");
				Test.printLog();
				System.exit(0);
			}
			File ffam = new File(fam);
			if (!ffam.exists()) {
				System.err.println("could not open " + fam +".");
				Test.LOG.append("could not open " + fam + ".\n");
				Test.printLog();
				System.exit(0);
			}
			bfileFlag = true;
		}

		if (cl.hasOption(cmd_covar)) {
			pheno = cl.getOptionValue(cmd_covar);
			File fpheno = new File(pheno);
			if (!fpheno.exists()) {
				System.err.println("could not open " + fpheno + ".");
				Test.LOG.append("could not open " + fpheno + ".\n");
				Test.printLog();
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_header)) {
			header = true;
		}
		// if (cl.hasOption(cmd_topN)) {
		// topN = Integer.parseInt(cl.getOptionValue(cmd_topN));
		// }
		if (cl.hasOption(cmd_pheno_number)) {
			response = Integer.parseInt(cl.getOptionValue(cmd_pheno_number)) - 1;
			if (response < -1) {
				System.err.println("bad parameter for --" + cmd_pheno_number_long + ": " +(response+1)+ ".");
				Test.LOG.append("bad parameter for --" + cmd_pheno_number_long +  ": " +(response+1)+ ".\n");
				Test.printLog();
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_covar_number)) {
			String[] p = cl.getOptionValues(cmd_covar_number);
			HashSet<Integer> idx = NewIt.newHashSet();
			for (int i = 0, len = p.length; i < len; i++) {
				if (p[i].contains("-")) {
					String[] pp = p[i].split("-");
					if (pp.length != 2) {
						System.err.println("bad parameter for option --" + cmd_covar_number_long + ": " + p[i] +".");
						Test.LOG.append("bad parameter for option --" + cmd_covar_number_long + ": " + p[i] +".\n");
						Test.printLog();
						System.exit(0);
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
					System.err.println("bad parameter for option --" + cmd_covar_number_long + ": " + predictor[c] +".");
					Test.LOG.append("bad parameter for option --" + cmd_covar_number_long + ": " + predictor[c] +".\n");
					Test.printLog();
					System.exit(0);
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
						System.err.println("bad parameter for option --" + cmd_covar_name_long + ": " + p[i] + ".");
						Test.LOG.append("bad parameter for option --" + cmd_covar_name_long + ": " + p[i] + ".\n");
						Test.printLog();
						System.exit(0);
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

		if (cl.hasOption(cmd_bgsnp)) {
			String[] bg = cl.getOptionValues(cmd_bgsnp);
			HashSet<String> bgSet = NewIt.newHashSet();
			for (int i = 0; i < bg.length; i++) {
				bgSet.add(bg[i]);
			}
			if (bgSet.size() != bg.length) {
				System.err.println("bad parameter for --" + cmd_bgsnp + ".");
				Test.LOG.append("bad parameter for --" + cmd_bgsnp + ".\n");
				Test.printLog();
				System.exit(0);
			}
			bgsnp = cl.getOptionValues(cmd_bgsnp);
			bgsnpFlag = true;
		}

		
		if (cl.hasOption(cmd_hg18)) {
			hg18Flag = true;
			hg19Flag = false;
			hgFile = "/gene36.txt";
		}

		if (cl.hasOption(cmd_hg19)) {
			hg19Flag = true;
			hg18Flag = false;
			hgFile = "/gene37.txt";
		}

		if (cl.hasOption(cmd_hg)) {
			hg19Flag = false;
			hg18Flag = false;
			hgFile = cl.getOptionValue(cmd_hg);
		}

		if (cl.hasOption(cmd_snp2genelist)) {
			snp2genefileFlag = true;
		}

		if (cl.hasOption(cmd_snp2genemlist)) {
			snp2genefilesFlag = true;
		}

		if (cl.hasOption(cmd_region)) {
			String[] r = cl.getOptionValues(cmd_region);
			ArrayList<String> chr = NewIt.newArrayList();
			ArrayList<String> b = NewIt.newArrayList();
			ArrayList<String> e = NewIt.newArrayList();

			for (int i = 0; i < r.length; i++) {
				String[] s = r[i].split(",");
				if (s.length != 3) {
					System.err.println("bad parameter for --" + cmd_region + ": " + r[i] +".");
					Test.LOG.append("bad parameter for --" + cmd_region + ": " + r[i] +".\n");
					Test.printLog();
					System.exit(0);
				}
				chr.add(s[0]);
				b.add(s[1]);
				e.add(s[2]);
			}
			chr_reg = (String[]) chr.toArray(new String[0]);
			begin = new double[b.size()];
			end = new double[e.size()];
			for (int i = 0; i < r.length; i++) {
				begin[i] = Double.parseDouble(b.get(i));
				end[i] = Double.parseDouble(e.get(i));
			}
			regionFlag = true;
		}

		if (cl.hasOption(cmd_gene_list)) {
			String gl = cl.getOptionValue(cmd_gene_list);

			File f = new File(gl);
			if (!f.exists()) {
				System.err.println("could not find file for --option " + cmd_gene_list_long + ": " + gl +".");
				Test.LOG.append("could not find file for --option " + cmd_gene_list_long + ": " + gl + ".\n");
				Test.printLog();
				System.exit(0);
			}
			BufferedReader reader0 = null;
			try {
				reader0 = new BufferedReader(new FileReader(f));
			} catch (IOException E) {
				System.err.println("could not open gene list " + gl + ".");
				Test.LOG.append("could not open gene list " + gl + ".\n");
				Test.printLog();
				System.exit(0);
			}
			String line0 = null;
			HashSet<String> gSet = NewIt.newHashSet();
			try {
				while((line0 = reader0.readLine()) != null) {
					String[] gn = line0.split(delim);
					gSet.add(gn[0]);
				}
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			String[] g = (String[]) gSet.toArray(new String[0]);
			boolean[] gflag = new boolean[gSet.size()];
			Arrays.fill(gflag, false);
			ArrayList<String> ge = NewIt.newArrayList();
			ArrayList<String> g_chr = NewIt.newArrayList();
			ArrayList<String> g_begin = NewIt.newArrayList();
			ArrayList<String> g_end = NewIt.newArrayList();
			
			BufferedReader reader = null;
			if(hg18Flag || hg19Flag) {
				InputStream is = getClass().getResourceAsStream(hgFile);
				DataInputStream in = new DataInputStream(is);
				reader = new BufferedReader(new InputStreamReader(in));
			} else {
				File fhg = new File(hgFile);
				if (!fhg.exists()) {
					System.err.println("could not find file for --option " + cmd_hg + ": " + hgFile +".");
					Test.LOG.append("could not find file for --option " + cmd_hg + ": " + hgFile + ".\n");
					Test.printLog();
					System.exit(0);
				}

				try {
					reader = new BufferedReader(new FileReader(fhg));
				} catch (IOException E) {
					System.err.println("could not open gene list " + hgFile + ".");
					Test.LOG.append("could not open gene list " + hgFile + ".\n");
					Test.printLog();
					System.exit(0);
				}

			}

			String line = null;
			try {
				while ((line = reader.readLine()) != null) {

					String[] s = line.split("\\s+");
//					System.err.println(line);
					if (s.length != 4) {
						continue;
					}

					for (int i = 0; i < g.length; i++) {
						if (s[0].compareTo(g[i]) == 0) {
							ge.add(s[0]);
							g_chr.add(s[1]);
							g_begin.add(s[2]);
							g_end.add(s[3]);
							gflag[i] = true;
						}
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			boolean flag = true;
			int count = 0;
			for(int i = 0; i < gflag.length; i++) {
				if(!gflag[i]) {
					System.err.println("could not fine gene " + g[i] + ".");
					Test.LOG.append("could not find gene " + g[i] + ".\n");
					flag = false;
					count ++;
				}
			}

			System.err.println("of " + gflag.length + " genes " + (gflag.length - count) + " was found.");
			Test.LOG.append("of " + gflag.length + " genes " + (gflag.length - count) + " was found.\n");
			
			if (!snp2genefileFlag && !snp2genefilesFlag) {
				if(!flag) {
					Test.printLog();
					System.exit(0);
				}
			}

			gene = (String[]) ge.toArray(new String[0]);
			gene_chr = (String[]) g_chr.toArray(new String[0]);
			gene_begin = new double[gene_chr.length];
			gene_end = new double[gene_chr.length];

			for (int i = 0; i < gene_chr.length; i++) {
				gene_begin[i] = Double.parseDouble(g_begin.get(i)) / 1000;
				gene_end[i] = Double.parseDouble(g_end.get(i)) / 1000;
				System.err.println(gene[i] + ": chr" + gene_chr[i] + " " +gene_begin[i] + "k ~ " + gene_end[i] + "k.");
				Test.LOG.append(gene[i] + ": chr" + gene_chr[i] + " " +gene_begin[i] + "k ~ " + gene_end[i] + "k.\n");
			}
			geneFlag = true;
		}
		
		if (cl.hasOption(cmd_gene_window)) {
			double gw = Double.parseDouble(cl.getOptionValue(cmd_gene_window));
			if (gw < 0) {
				System.err.println("bad parameter for option --" + cmd_gene_window_long + ": " + gw + ".");
				Test.LOG.append("bad parameter for option --" + cmd_gene_window_long + ": " + gw +".\n");
				Test.printLog();
				System.exit(0);
			}
			genewindow = gw;
		}

		if (cl.hasOption(cmd_gene)) {
			String[] g = cl.getOptionValues(cmd_gene);
			boolean[] gflag = new boolean[g.length];
			Arrays.fill(gflag, false);
			ArrayList<String> ge = NewIt.newArrayList();
			ArrayList<String> g_chr = NewIt.newArrayList();
			ArrayList<String> g_begin = NewIt.newArrayList();
			ArrayList<String> g_end = NewIt.newArrayList();
			
			BufferedReader reader = null;
			if(hg18Flag || hg19Flag) {
				InputStream is = getClass().getResourceAsStream(hgFile);
				DataInputStream in = new DataInputStream(is);
				reader = new BufferedReader(new InputStreamReader(in));
			} else {
				File fhg = new File(hgFile);
				if (!fhg.exists()) {
					System.err.println("could not find file for --option " + cmd_hg + ": " + hgFile +".");
					Test.LOG.append("could not find file for --option " + cmd_hg + ": " + hgFile + ".\n");
					Test.printLog();
					System.exit(0);
				}
				try {
					reader = new BufferedReader(new FileReader(fhg));
				} catch (IOException E) {
					System.err.println("could not open gene list " + hgFile + ".");
					Test.LOG.append("could not open gene list " + hgFile + ".\n");
					Test.printLog();
					System.exit(0);
				}

			}

			String line = null;
			try {
				while ((line = reader.readLine()) != null) {

					String[] s = line.split("\\s+");
//					System.err.println(line);
					if (s.length != 4) {
						continue;
					}

					for (int i = 0; i < g.length; i++) {
						if (s[0].compareTo(g[i]) == 0) {
							ge.add(s[0]);
							g_chr.add(s[1]);
							g_begin.add(s[2]);
							g_end.add(s[3]);
							gflag[i] = true;
						}
					}

				}
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			boolean flag = true;
			int count = 0;
			for(int i = 0; i < gflag.length; i++) {
				if(!gflag[i]) {
					System.err.println("could not find gene " + g[i] + ".");
					Test.LOG.append("could not find gene " + g[i] + ".\n");
					flag = false;
					count++;
				}
			}
			System.err.println("of " + gflag.length + " genes " + (gflag.length - count) + " was found.");
			Test.LOG.append("of " + gflag.length + " genes " + (gflag.length - count) + " was found.\n");
			if (!snp2genefileFlag && !snp2genefilesFlag) {
				if(!flag) {
					Test.printLog();
					System.exit(0);
				}
			}

			gene = (String[]) ge.toArray(new String[0]);
			gene_chr = (String[]) g_chr.toArray(new String[0]);
			gene_begin = new double[gene_chr.length];
			gene_end = new double[gene_chr.length];

			for (int i = 0; i < gene_chr.length; i++) {
				gene_begin[i] = Double.parseDouble(g_begin.get(i)) / 1000;
				gene_end[i] = Double.parseDouble(g_end.get(i)) / 1000;
				System.err.println(gene[i] + ": chr" + gene_chr[i] + " " +gene_begin[i] + "k ~ " + gene_end[i] + "k");
				Test.LOG.append(gene[i] + ": chr" + gene_chr[i] + " " +gene_begin[i] + "k ~ " + gene_end[i] + "k.\n");
			}
			geneFlag = true;
		}

		if (cl.hasOption(cmd_extract)) {

			if (!transFlag) {
				String[] snps_file = cl.getOptionValues(cmd_extract);
				ArrayList<String> includesnpList = NewIt.newArrayList();
				for (int h = 0; h < snps_file.length; h++) {
					File f = new File(snps_file[h]);
					if (!f.exists()) {
						System.err.println("could not find --" + cmd_extract + ": " + snps_file[h] + ".");
						Test.LOG.append("could not fine --" + cmd_extract + ": " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					BufferedReader reader = null;
					try {
						reader = new BufferedReader(new FileReader(f));
					} catch (IOException E) {
						System.err.println("could not read --" + cmd_extract + ": " + snps_file[h] + ".");
						Test.LOG.append("could not read --" + cmd_extract + ": " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					ArrayList<String> snp = NewIt.newArrayList();
					String line = null;
					try {
						while ((line = reader.readLine()) != null) {
							String[] s = line.split(delim);
							snp.add(s[0]);
						}
						reader.close();
					} catch (IOException E) {
						System.err.println("bad lines in " + snps_file[h] + ".");
						Test.LOG.append("bad lines in " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					if (snp.size() > 0) {
						ArrayList<String> insnp = NewIt.newArrayList();
						for (int i = 0; i < snp.size(); i++) {
							String subSNP = snp.get(i);
							insnp.add(subSNP);
						}

						if (insnp.size() > 0) {
							includesnpList.addAll(insnp);
							snpFlag = true;
						}
					}
					if (includesnpList.size() > 0) {
						includesnp = (String[]) includesnpList
								.toArray(new String[0]);
					}

				}
			} else {
				String[] snps_file = cl.getOptionValues(cmd_extract);
				xincludesnp = new String[snps_file.length][];
				for (int h = 0; h < snps_file.length; h++) {
					File f = new File(snps_file[h]);
					if (!f.exists()) {
						System.err.println("could not find " + snps_file[h] + ".");
						Test.LOG.append("could not find " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					BufferedReader reader = null;
					try {
						reader = new BufferedReader(new FileReader(f));
					} catch (IOException E) {

						System.err.println("could not read " + snps_file[h] + ".");
						Test.LOG.append("could not read " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					ArrayList<String> snp = NewIt.newArrayList();
					String line = null;
					try {
						while ((line = reader.readLine()) != null) {
							String[] s = line.split(delim);
							snp.add(s[0]);
						}
						reader.close();
					} catch (IOException E) {

						System.err.println("bad lines in " + snps_file[h] + ".");
						Test.LOG.append("bad lines in " + snps_file[h] + ".\n");
						Test.printLog();
						System.exit(0);
					}
					if (snp.size() > 0) {
						ArrayList<String> insnp = NewIt.newArrayList();

						for (int i = 0; i < snp.size(); i++) {
							String subSNP = snp.get(i);
							insnp.add(subSNP);
						}
						if (insnp.size() > 0) {
							xincludesnp[h] = (String[]) insnp
									.toArray(new String[0]);
							snpFlag = true;
						}
					}
				}
			}
		}


		if (cl.hasOption(cmd_exclude)) {

			if (!transFlag) {
				String snps_file = cl.getOptionValue(cmd_exclude);
				ArrayList<String> excludesnpList = NewIt.newArrayList();
				for (int h = 0; h < 1; h++) {
					File f = new File(snps_file);
					if (!f.exists()) {
						System.err.println("could not find --" + cmd_extract + ": " + snps_file + ".");
						Test.LOG.append("could not fine --" + cmd_extract + ": " + snps_file + ".\n");
						Test.printLog();
						System.exit(0);
					}
					BufferedReader reader = null;
					try {
						reader = new BufferedReader(new FileReader(f));
					} catch (IOException E) {
						System.err.println("could not read --" + cmd_extract + ": " + snps_file + ".");
						Test.LOG.append("could not read --" + cmd_extract + ": " + snps_file + ".\n");
						Test.printLog();
						System.exit(0);
					}
					ArrayList<String> snp = NewIt.newArrayList();
					String line = null;
					try {
						while ((line = reader.readLine()) != null) {
							String[] s = line.split(delim);
							snp.add(s[0]);
						}
						reader.close();
					} catch (IOException E) {
						System.err.println("bad lines in " + snps_file + ".");
						Test.LOG.append("bad lines in " + snps_file + ".\n");
						Test.printLog();
						System.exit(0);
					}
					if (snp.size() > 0) {
						ArrayList<String> exsnp = NewIt.newArrayList();
						for (int i = 0; i < snp.size(); i++) {
							String subSNP = snp.get(i);
							exsnp.add(subSNP);
						}
						if (exsnp.size() > 0) {
							excludesnpList.addAll(exsnp);
							snpFlag = true;
						}
					}
					if (excludesnpList.size() > 0) {
						excludesnp = (String[]) excludesnpList
								.toArray(new String[0]);
					}
				}
			}
			
		}

		
		if (cl.hasOption(cmd_keep_male)) {
			keep_maleFlag = true;
		}
		if (cl.hasOption(cmd_keep_female)) {
			keep_femaleFlag = true;
		}
		if (cl.hasOption(cmd_ex_nosex)) {
			ex_nosexFlag = true;
		}
		if (cl.hasOption(cmd_remove)) {
			String file = cl.getOptionValue(cmd_remove);
			File f = new File(file);
			if (!f.exists()) {
				System.err.println("could not open " + file + ".");
				Test.LOG.append("could not open " + file + ".\n");
				Test.printLog();
				System.exit(0);
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(new File(file)));
			} catch (IOException E) {
				System.err.println("could not read " + file + ".");
				Test.LOG.append("coudl not read " + file + ".\n");
				Test.printLog();
				System.exit(0);
			}
			ArrayList<String> famList = NewIt.newArrayList();
			ArrayList<String> indList = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					String[] l = line.split(MDRConstant.delim);
					if(l.length < 2) continue;
					famList.add(l[0]);
					indList.add(l[1]);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			ex_family = new String[2][];
			ex_family[0] = (String[]) famList.toArray(new String[0]);
			ex_family[1] = (String[]) indList.toArray(new String[0]);
			removeFlag = true;
		}

		if (cl.hasOption(cmd_keep)) {
			String file = cl.getOptionValue(cmd_keep);
			File f = new File(file);
			if (!f.exists()) {
				System.err.println("could not open " + file + ".");
				Test.LOG.append("could not open " + file + ".\n");
				Test.printLog();
				System.exit(0);
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(new File(file)));
			} catch (IOException E) {
				System.err.println("could not read " + file + ".");
				Test.LOG.append("coudl not read " + file + ".\n");
				Test.printLog();
				System.exit(0);
			}
			ArrayList<String> famList = NewIt.newArrayList();
			ArrayList<String> indList = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					String[] l = line.split(MDRConstant.delim);
					if(l.length < 2) continue;
					famList.add(l[0]);
					indList.add(l[1]);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			indKeep = new String[2][];
			indKeep[0] = (String[]) famList.toArray(new String[0]);
			indKeep[1] = (String[]) indList.toArray(new String[0]);
			keepFlag = true;
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
				System.err.println("bad parameter for optin --" + cmd_chr + ".");
				Test.LOG.append("bad parameter for option --" + cmd_chr + ".\n");
				Test.printLog();
				System.exit(0);
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
					System.err.println("bad parameter for optin --" + cmd_snpwindow_long + " " + s[i] + ".");
					Test.LOG.append("bad parameter for option --" + cmd_snpwindow_long + " " + s[i] + ".\n");
					Test.printLog();
					System.exit(0);
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

				System.err.println("bad parameter for optin --" + cmd_maf + " " + maf + ".");
				Test.LOG.append("bad parameter for option --" + cmd_maf + " " + maf + ".\n");
				Test.printLog();
				System.exit(0);
			}
			mafFlag = true;
		}

		if (cl.hasOption(cmd_geno)) {
			geno = Double.parseDouble(cl.getOptionValue(cmd_geno));
			if (geno < 0) {

				System.err.println("bad parameter for optin --" + cmd_geno + " " + geno + ".");
				Test.LOG.append("bad parameter for option --" + cmd_geno + " " + geno + ".\n");
				Test.printLog();
				System.exit(0);
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

		if (cl.hasOption(cmd_reg)) {
			linkfunction = Integer.parseInt(cl.getOptionValue(cmd_reg));
		}
		if (cl.hasOption(cmd_cv)) {
			cv = Integer.parseInt(cl.getOptionValue(cmd_cv));
			if (cv < 2) {

				System.err.println("bad parameter for optin --" + cmd_cv + " " + cv + ".");
				Test.LOG.append("bad parameter for option --" + cmd_cv + " " + cv + ".\n");
				Test.printLog();
				System.exit(0);
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

				System.err.println("bad parameter for optin --" + ep + " " + ep + ".");
				Test.LOG.append("bad parameter for option --" + ep + " " + ep + ".\n");
				Test.printLog();
				System.exit(0);
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
		}
		if (cl.hasOption(cmd_thin)) {
			thin = Double.parseDouble(cl.getOptionValue(cmd_thin));
			if (thin < 0) {

				System.err.println("bad parameter for optin --" + cmd_thin + " " + thin + ".");
				Test.LOG.append("bad parameter for option --" + cmd_thin + " " + thin + ".\n");
				Test.printLog();
				System.exit(0);
			}
		}
		if (cl.hasOption(cmd_slice)) {
			String[] s = cl.getOptionValue(cmd_slice).split("/");
			slice = Integer.parseInt(s[0]);
			sliceN = Integer.parseInt(s[1]);
			if (slice <= 0 || sliceN <= 0 || slice > sliceN) {

				System.err.println("bad parameter for optin --" + cmd_slice + " " + slice + ".");
				Test.LOG.append("bad parameter for option --" + cmd_slice + " " + slice + ".\n");
				Test.printLog();
				System.exit(0);
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
			status_shift = 0;
			status_shiftFlag = true;
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

		if (cl.hasOption(cmd_verbose)) {
			verboseFlag = true;
		}
		if (cl.hasOption(cmd_version)) {
			System.err.println();
			Test.printLog();
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
