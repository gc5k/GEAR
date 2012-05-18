package parameter;

import java.io.File;
import java.util.HashSet;
import java.util.Random;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import test.Test;
import util.NewIt;

public class Parameter {

	private final String sep=",";
	public static String whitespace = "\\s+";

////////////////bfile
	private final String cmd_bfile = "bfile";
	public static String bfile = null;
	public static String bedfile = null;
	public static String bimfile = null;
	public static String famfile = null;
	public static boolean bfileOption = false;

	private final String cmd_file = "file";
	public static String tbfile = null;
	public static String pedfile = null;
	public static String mapfile = null;
	public static boolean fileOption = false;

	public static boolean header = false;

///////////////real check
	private final String cmd_realcheck_threshold_upper = "realcheck_threshold_upper";
	private final String cmd_realcheck_threshold_upper_long = "realcheck-threshold-upper";
	public static double realcheckThresholdUpper = 1;
	
	private final String cmd_realcheck_threshold_lower = "realcheck_threshold_lower";
	private final String cmd_realcheck_threshold_lower_long = "realcheck-threshold-lower";
	public static double realcheckThresholdLower = 0;	
	public static boolean realcheckThresholdUpperLowerFlag = false;

	private final String cmd_realcheck_threshold = "realcheck_threshold";
	private final String cmd_realcheck_threshold_long = "realcheck-threshold";
	public static double realcheckThreshold = 1;
	public static boolean realcheckThresholdFlag = true;

	private final String cmd_realcheck_marker = "realcheck_marker";
	private final String cmd_realcheck_marker_long = "realcheck-marker-number";
	public static int realcheckMarkerNumber = 0;

	private final String cmd_realcheck = "realcheck";
	public static boolean realcheckFlag = false;

	private final String cmd_bfile2 = "bfile2";
	public static String bfile2 = null;
	public static String bedfile2 = null;
	public static String bimfile2 = null;
	public static String famfile2 = null;
	public static boolean bfileOption2 = false;

///////////////simulation nuclear family
	private final String cmd_simu_fam = "simu_fam";
	private final String cmd_simu_fam_long = "simu-fam";
	public static boolean simufamFlag = false;
	
	private final String cmd_simu_fam_size = "simu_fam_size";
	private final String cmd_simu_fam_size_long = "simu-fam-size";
	public static int simu_fam_size = 100;
	private final String cmd_simu_fam_marker = "simu_fam_marker";
	private final String cmd_simu_fam_marker_long = "simu-fam-marker";
	public static int simu_fam_marker = 10;

///////////////simulation real data
	private final String cmd_simu_seed = "simu_seed";
	private final String cmd_simu_seed_long = "simu-seed";
	public static long simuSeed = (new Random()).nextLong();

	private final String cmd_simu_rep = "simu_rep";
	private final String cmd_simu_rep_long = "simu-rep";
	public static int simuRep = 1;

	private final String cmd_simu_casual_loci = "simu_casual_loci";
	private final String cmd_simu_casual_loci_long = "simu-casual-loci";
	public static String simuCasualLoci = null;

	private final String cmd_simu_rnd_casual_loci = "simu_rnd_casual_loci";
	private final String cmd_simu_rnd_casual_loci_long = "simu-rnd-casual-loci";
	public static int simuRndCasualLoci = 0;

	private final String cmd_simu_hsq = "simu_hsq";
	private final String cmd_simu_hsq_long = "simu-hsq";
	public static double simuHsq = 0.5;

	private final String cmd_simu_qt = "simu_qt";
	private final String cmd_simu_qt_long = "simu-qt";

	private final String cmd_simu_cc = "simu_cc";
	private final String cmd_simu_cc_long = "simu-cc";
	public static int[] simuCC = { 0, 0 };

	private final String cmd_simu_k = "simu_k";
	private final String cmd_simu_k_long = "simu-k";
	public static double simuK = 0.1;

	public final static int sm_qt = 0;
	public final static int sm_cc = 1;
	public static boolean simuFlag = false;
	public static boolean[] simuType = { true, false };

/////////////////simulation polygenic
	private final String cmd_simu_poly_qt = "simu_poly_qt";
	private final String cmd_simu_poly_qt_long = "simu-poly-qt";	
	public static boolean simupolyQTFlag = false;

	private final String cmd_simu_poly_cc = "simu_poly_cc";
	private final String cmd_simu_poly_cc_long = "simu-poly-cc";

	public static boolean simupolyCCFlag = false;

	private final String cmd_poly_loci = "poly_loci";
	private final String cmd_poly_loci_long = "poly-loci";
	
	public static int polyLoci = 1000;
	
	private final String cmd_poly_LD = "poly_ld";
	private final String cmd_poly_LD_long = "poly-ld";
	public static double polyLD = 0;

	private final String cmd_poly_U = "poly_U";
	private final String cmd_poly_U_long = "poly-U";
	public static boolean polyU = false;

	private final String cmd_poly_sample = "poly_sample_size";
	private final String cmd_poly_sample_long = "poly-sample-size";
	public static int polySample = 1000;

	private final String cmd_poly_case = "poly_case";
	private final String cmd_poly_case_long = "poly-case";	
	public static int polyCS = 500;

	private final String cmd_poly_k = "poly_k";
	private final String cmd_poly_k_long = "poly-k";
	public static double polyK = 0.05;

	private final String cmd_poly_hsq = "poly_hsq";
	private final String cmd_poly_hsq_long = "poly-hsq";
	public static double polyHsq = 0.5;

///////////////////nontrans
	private final String cmd_nontrans = "nontrans";
	public static boolean nontransFlag = false;

	private final String cmd_nontrans_seed = "nontrans_seed";
	private final String cmd_nontrans_seed_long = "nontrans-seed";	
	public static long nontransSeed = 2010;

	private final String cmd_nontrans_cases = "nontrans_cases";
	private final String cmd_nontrans_cases_long = "nontrans-cases";
	public static boolean nontranscasesFlag = false;

	private final String cmd_nontrans_controls = "nontrans_controls";
	private final String cmd_nontrans_controls_long = "nontrans-controls";
	public static boolean nontranscontrolsFlag = false;

///////////////////pop stat
	private final String cmd_freq = "freq";
	private final String cmd_geno_freq = "geno_freq";
	private final String cmd_geno_freq_long = "geno-freq";

	public final String cmd_sum_stat_help = "sum_stat_help";
	public final String cmd_sum_stat_help_long = "sum-stat-help";
	public final static int freq = 0;
	public final static int geno_freq = 1;
	public static boolean sumStatFlag = false;
	public static boolean[] sumStatOption = { true, false };

//he regression
	public static boolean heFlag = true;
	private final String cmd_he = "he";
	private final String cmd_grm = "grm";
	public static String grm = null;
	public static String grm_id = null;

	private final String cmd_pheno = "pheno";
	public static String pheno = null;

	private final String cmd_mpheno = "mpheno";
	public static int[] mpheno = {1};

	private final String cmd_reverse = "reverse";
	public static boolean reverse = false;

	private final String cmd_k = "k";
	public static boolean k_button = false;
	public static double k = 0.01;
	
	private final String cmd_he_sd = "he_sd";
	private final String cmd_he_sd_long = "he-sd";
	private final String cmd_he_ss = "he_ss";
	private final String cmd_he_ss_long = "he-ss";
	private final String cmd_he_cp = "he_cp";
	private final String cmd_he_cp_long = "he-cp";
	public static int he_sd = 0;
	public static int he_ss = 1;
	public static int he_cp = 2;
	public static boolean[] heType = {true, false, false};

	private final String cmd_out = "out";
	public static String out = "he";

///////////////////heritability transformation
	public static boolean calOption = false;
	private final String cmd_cal_k = "cal_k";
	private final String cmd_cal_k_long = "cal-k";
	public double cal_k = 0;
	
	private final String cmd_cal_hl = "cal_hl";
	private final String cmd_cal_hl_long = "cal-hl";
	public double cal_hl = 0;
	
	private final String cmd_cal_ho = "cal_ho";
	private final String cmd_cal_ho_long = "cal-ho";
	public double cal_ho = 0;
	
	private final String cmd_cal_h2_se = "cal_h2_se";
	private final String cmd_cal_h2_se_long = "cal-h2-se";
	public double cal_h2_se = -1;
	
	private final String cmd_cal_cc = "cal_cc";
	private final String cmd_cal_cc_long = "cal-cc";
	public double[] cal_cc = {0,0}; 

	private final String cmd_na = "na";
	public static String[] na = {"-9", "-NA"};

///////////////////level 1 snp selection
	private final String cmd_chr = "chr";
	public static String[] inchr = null;
	public static String[] exchr = null;
	public static boolean inchrFlag = false;
	public static boolean exchrFlag = false;

//////////////////level 1 individual selection
	private final String cmd_keep = "keep";
	public static String keepFile = null;
	public static boolean keepFlag = false;

	// individual selection start
	private final String cmd_remove = "remove";
	public static String removeFile = null;
	public static boolean removeFlag = false;

/////////////////write bed file
	private final String cmd_make_bed = "make_bed";
	public static String cmd_make_bed_long = "make-bed";
	public static boolean makebedFlag = false;

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
///////////////////global	

	
	public static String version = "0.1, April/06/2012\n";

	public static boolean transFlag = false;

	public static String missing_allele = "0";

	public static boolean status_shiftFlag = false;

	public static String missing_phenotype = "-9";

	public static double status_shift = -1;

	public static String missingGenotype = "22";

	public static boolean covar_header_flag = false;

	public static long seed = 2012;

	public static boolean genoFlag = false;

	public static double geno = 0;

	public static boolean mafFlag = false;
	
	public static double maf = 0;

	public static boolean maxmafFlag = false;
	
	public static double max_maf = 0.55;

	private final String cmd_help = "help";

	
	private Options ops = new Options();

	private CommandLineParser parser = new PosixParser();

	public Parameter() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	public void commandInitial() {

		ops.addOption(OptionBuilder.withDescription("file ").hasArg().create(cmd_file));

		ops.addOption(OptionBuilder.withDescription("bfile ").hasArg().create(cmd_bfile));

//realcheck
		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_upper_long).withDescription("realcheck marker threshold upper bounder").hasArg().create(cmd_realcheck_threshold_upper));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_lower_long).withDescription("realcheck marker threshold lower bounder").hasArg().create(cmd_realcheck_threshold_lower));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_long).withDescription("realcheck marker threshold").hasArg().create(cmd_realcheck_threshold));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_marker_long).withDescription("realcheck marker number").hasArg().create(cmd_realcheck_marker));

		ops.addOption(OptionBuilder.withDescription("realcheck ").create(cmd_realcheck));

		ops.addOption(OptionBuilder.withDescription("bfile2 ").hasArg().create(cmd_bfile2));

//simulation nuclear fam
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_long).withDescription("simulation nuclear family ").create(cmd_simu_fam));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_size_long).withDescription("simulation nuclear family size ").hasArg().create(cmd_simu_fam_size));
		
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_marker_long).withDescription("simulation number for nuclear family ").hasArg().create(cmd_simu_fam_marker));		
		
//simulation real data
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_seed_long).withDescription("gwas simulation seed ").hasArg().create(cmd_simu_seed));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_rep_long).withDescription("gwas simulation replication ").hasArg().create(cmd_simu_rep));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_casual_loci_long).withDescription("gwas simulation casual loci ").hasArg().create(cmd_simu_casual_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_rnd_casual_loci_long).withDescription("gwas simulation casual loci number ").hasArg().create(cmd_simu_rnd_casual_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_hsq_long).withDescription("gwas simulation heritability ").hasArg().create(cmd_simu_hsq));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_qt_long).withDescription("gwas simulate quantitative traits ").create(cmd_simu_qt));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_cc_long).withDescription("gwas simulate case-control ").hasArgs(2).create(cmd_simu_cc));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_k_long).withDescription("gwas prevalence of the binary trait ").hasArg().create(cmd_simu_k));

//nontransmitted
		ops.addOption(OptionBuilder.withDescription("nontransmitted ").create(cmd_nontrans));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_cases_long).withDescription("nontransmitted filter cases ").create(cmd_nontrans_cases));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_controls_long).withDescription("nontransmitted filter controls ").create(cmd_nontrans_controls));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_seed_long).withDescription("gwas prevalence of the binary trait ").hasArg().create(cmd_nontrans_seed));
		
//simulation polygenic model
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_poly_qt_long).withDescription("polygenic model simulation for quantitative traits ").create(cmd_simu_poly_qt));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_poly_cc_long).withDescription("polygenic model simulation for case-control sample ").create(cmd_simu_poly_cc));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_loci_long).withDescription("number of polygenic loci, defualt= " + polyLoci).hasArg().create(cmd_poly_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_LD_long).withDescription("LD (correlation), defualt= " + polyLD).hasArg().create(cmd_poly_LD));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_U_long).withDescription("polygenic model has Uniform Effect? " + polyU).create(cmd_poly_U));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_sample_long).withDescription("polygenic model sample size, default = " + polySample).hasArg().create(cmd_poly_sample));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_case_long).withDescription("number of cases, default = " + polyCS).hasArg().create(cmd_poly_case));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_k_long).withDescription("prevalence, default = " + polyK).hasArg().create(cmd_poly_k));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_hsq_long).withDescription("heritability, default = " + polyHsq).hasArg().create(cmd_poly_hsq));
		
//pop stat
		ops.addOption(OptionBuilder.withDescription("calculate MAF frequency ").create(cmd_freq));

		ops.addOption(OptionBuilder.withLongOpt(cmd_geno_freq_long).withDescription("calculate genotype frequency ").create(cmd_geno_freq));

//snp selection
		ops.addOption(OptionBuilder.withDescription("select chromosomes").hasArgs().create(cmd_chr));
//individual selection

		ops.addOption(OptionBuilder.withDescription("remove individuals").hasArg().create(cmd_remove));

		ops.addOption(OptionBuilder.withDescription("keep individuals").hasArg().create(cmd_keep));

		ops.addOption(OptionBuilder.withDescription("keep males only").withLongOpt(cmd_keep_male_long)
				.create(cmd_keep_male));
		ops.addOption(OptionBuilder.withDescription("keep females only").withLongOpt(cmd_keep_female_long)
				.create(cmd_keep_female));
		ops.addOption(OptionBuilder.withDescription("exclude unknown sex").withLongOpt(cmd_ex_nosex_long)
				.create(cmd_ex_nosex));

//make bed
		
		ops.addOption(OptionBuilder.withLongOpt(cmd_make_bed_long).withDescription("make bed ").create(cmd_make_bed));		
		
//haseman-elston regression
		ops.addOption(OptionBuilder.withDescription("haseman-elston regression ").hasArg().create(cmd_he));		
		
		ops.addOption(OptionBuilder.withDescription("grm ").hasArg().create(cmd_grm));

		ops.addOption(OptionBuilder.withDescription("pheno ").hasArg().create(cmd_pheno));

		ops.addOption(OptionBuilder.withLongOpt(cmd_mpheno).withDescription("pheno number " + cmd_mpheno).hasArg().withArgName("index").create(cmd_mpheno));

		ops.addOption(OptionBuilder.withDescription("reverse ").create(cmd_reverse));

		ops.addOption(OptionBuilder.withLongOpt(cmd_he_sd_long).withDescription("phenotype is coded as squared difference").create(cmd_he_sd));

		ops.addOption(OptionBuilder.withLongOpt(cmd_he_ss_long).withDescription("phenotype is coded as squared sum").create(cmd_he_ss));

		ops.addOption(OptionBuilder.withLongOpt(cmd_he_cp_long).withDescription("phenotype is coded as cross product").create(cmd_he_cp));

		ops.addOption(OptionBuilder.withDescription("prevalence ").hasArg().create(cmd_k));

		ops.addOption(OptionBuilder.withDescription("na ").hasArg().create(cmd_na));

		ops.addOption(OptionBuilder.withLongOpt(cmd_cal_k_long).withDescription("calculate heritability on the liability/observed scale with value K " + cmd_cal_k).hasArg().create(cmd_cal_k));
		
		ops.addOption(OptionBuilder.withLongOpt(cmd_cal_hl_long).withDescription("calculate heritability on the liability/observed scale " + cmd_cal_hl_long).hasArg().create(cmd_cal_hl));

		ops.addOption(OptionBuilder.withLongOpt(cmd_cal_ho_long).withDescription("calculate heritability on the liability/observed scale " + cmd_cal_ho_long).hasArg().create(cmd_cal_ho));

		ops.addOption(OptionBuilder.withLongOpt(cmd_cal_cc_long).withDescription("number of case and controls " + cmd_cal_cc).hasArg().create(cmd_cal_cc));

		ops.addOption(OptionBuilder.withLongOpt(cmd_cal_h2_se_long).withDescription("se of heritability on the liability/observed scale ").hasArg().create(cmd_cal_h2_se));

		ops.addOption(OptionBuilder.withDescription("root file, default = " + out).hasArg().create(cmd_out));

		ops.addOption(OptionBuilder.withDescription("help manual.").create(cmd_help));

	}

	public void commandListenor(String[] args) {
		CommandLine cl = null;
		try {
			cl = parser.parse(ops, args);
		} catch (ParseException E) {
			E.printStackTrace(System.err);
		}

		if (cl.hasOption(cmd_help)) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("HE Regression", ops);
			System.exit(1);
		}

		if (cl.hasOption(cmd_file)) {
			tbfile = cl.getOptionValue(cmd_file);
			pedfile = (new StringBuffer(tbfile)).append(".ped").toString();
			mapfile = (new StringBuffer(tbfile)).append(".map").toString();
			fileOption = true;
		}

		if (cl.hasOption(cmd_bfile)) {
			bfile = cl.getOptionValue(cmd_bfile);
			bedfile = (new StringBuffer(bfile)).append(".bed").toString();
			bimfile = (new StringBuffer(bfile)).append(".bim").toString();
			famfile = (new StringBuffer(bfile)).append(".fam").toString();
			bfileOption = true;
		}

//realcheck
		if (cl.hasOption(cmd_realcheck_threshold_upper)) {
			realcheckThresholdUpper = Double.parseDouble(cl.getOptionValue(cmd_realcheck_threshold_upper));
			if (realcheckThresholdUpper < 0 && realcheckThresholdUpper > 1) {
				System.err.println("realcheck threshold upper bounder should be tween 0 and 1");
				System.exit(0);
			}
			realcheckThresholdUpperLowerFlag = true;
		}
		if (cl.hasOption(cmd_realcheck_threshold_lower)) {
			realcheckThresholdLower = Double.parseDouble(cl.getOptionValue(cmd_realcheck_threshold_lower));
			if (realcheckThresholdLower < 0 && realcheckThresholdLower > 1) {
				System.err.println("realcheck threshold Lower bounder should be tween 0 and 1");
				System.exit(0);
			}
			realcheckThresholdUpperLowerFlag = true;
		}

		if (cl.hasOption(cmd_realcheck_threshold)) {
			realcheckThreshold = Double.parseDouble(cl.getOptionValue(cmd_realcheck_threshold));

			if (realcheckThreshold < 0) {
				System.err.println("realcheck threshold should be tween 0 and 1");
				System.exit(0);
			}
			realcheckThresholdFlag = true;
		}

		if (cl.hasOption(cmd_realcheck_threshold)) {
			realcheckThreshold = Double.parseDouble(cl.getOptionValue(cmd_realcheck_threshold));

			if (realcheckThreshold < 0) {
				System.err.println("realcheck threshold should be tween 0 and 1");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_realcheck_marker)) {
			realcheckMarkerNumber = Integer.parseInt(cl.getOptionValue(cmd_realcheck_marker));
			if (realcheckMarkerNumber < 0) {
				System.err.println("realcheck marker number should be greater than 0");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_realcheck)) {
			realcheckFlag = true;
		}

		if (cl.hasOption(cmd_bfile2)) {
			bfile2 = cl.getOptionValue(cmd_bfile2);
			bedfile2 = (new StringBuffer(bfile2)).append(".bed").toString();
			bimfile2 = (new StringBuffer(bfile2)).append(".bim").toString();
			famfile2 = (new StringBuffer(bfile2)).append(".fam").toString();
			bfileOption2 = true;
		}

//snp selection
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
				chr = (String[]) chrSet.toArray(new String[0]);
				inchrFlag = true;
			}
			if (exSet.size() > 0) {
				exchr = (String[]) exSet.toArray(new String[0]);
				exchrFlag = true;
			}
		}
		
//individual selection 1 keep
		if (cl.hasOption(cmd_keep)) {
			keepFile = cl.getOptionValue(cmd_keep);
			File f = new File(keepFile);
			if (!f.exists()) {
				System.err.println("could not open " + keepFile + ".");
				Test.LOG.append("could not open " + keepFile + ".\n");
				Test.printLog();
				System.exit(0);
			}
			keepFlag = true;
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
			removeFile = cl.getOptionValue(cmd_remove);
			File f = new File(removeFile);
			if (!f.exists()) {
				System.err.println("could not open " + removeFile + ".");
				Test.LOG.append("could not open " + removeFile + ".\n");
				Test.printLog();
				System.exit(0);
			}
			removeFlag = true;
		}

//make bed
		if (cl.hasOption(cmd_make_bed)) {
			makebedFlag = true;
		}

//nontrans		
		if (cl.hasOption(cmd_nontrans) ) {
			nontransFlag = true;
		}

		if (cl.hasOption(cmd_nontrans_seed)) {
			nontransSeed = Long.parseLong(cl.getOptionValue(cmd_nontrans_seed));
		}

		if (cl.hasOption(cmd_nontrans_cases)) {
			nontranscasesFlag = true;
			nontranscontrolsFlag = false;
		}

		if (cl.hasOption(cmd_nontrans_controls)) {
			nontranscontrolsFlag = true;
			nontranscasesFlag = false;
		}

//pop stat
		if (cl.hasOption(cmd_freq) ) {
			sumStatFlag = true;
			sumStatOption[freq] = true;
			sumStatOption[geno_freq] = false;
		}
		if (cl.hasOption(cmd_geno_freq)) {
			sumStatFlag = true;
			sumStatOption[freq] = false;
			sumStatOption[geno_freq] = true;
		}

//simulation nuclear family		

		if (cl.hasOption(cmd_simu_fam)) {
			simufamFlag = true;
		}
		
		if (cl.hasOption(cmd_simu_fam_size)) {
			simu_fam_size = Integer.parseInt(cl.getOptionValue(cmd_simu_fam_size));
		}
		
		if (cl.hasOption(cmd_simu_fam_marker)) {
			simu_fam_marker = Integer.parseInt(cl.getOptionValue(cmd_simu_fam_marker));
		}
		
//simulation real data
		if (cl.hasOption(cmd_simu_qt)) {
			simuFlag = true;
			simuType[sm_qt] = false;
			simuType[sm_cc] = true;
		}

		if (cl.hasOption(cmd_simu_cc)) {
			simuFlag = true;
			simuType[sm_qt] = true;
			simuType[sm_cc] = false;
			
			String[] s = cl.getOptionValues(cmd_simu_cc);
			simuCC[0] = Integer.parseInt(s[0]);
			simuCC[1] = Integer.parseInt(s[1]);
		}

		if (cl.hasOption(cmd_simu_seed)) {
			simuSeed = Long.parseLong(cl.getOptionValue(cmd_simu_seed));
		}

		if (cl.hasOption(cmd_simu_hsq)) {
			simuHsq = Double.parseDouble(cl.getOptionValue(cmd_simu_hsq));
			if (simuHsq < 0 || simuHsq > 1) {
				System.err
						.println("simulation heritability should be between 0 and 1");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_simu_casual_loci)) {
			simuCasualLoci = cl.getOptionValue(cmd_simu_casual_loci);
		}

		if (cl.hasOption(cmd_simu_rnd_casual_loci)) {
			simuRndCasualLoci = Integer.parseInt(cl.getOptionValue(cmd_simu_rnd_casual_loci));
		}

		if (cl.hasOption(cmd_simu_k)) {
			simuK = Double.parseDouble(cl.getOptionValue(cmd_simu_k));
			if (simuK < 0 || simuK > 1) {
				System.err
						.println("simulation prevalence should be between 0 and 1");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_simu_rep)) {
			simuRep = Integer.parseInt(cl.getOptionValue(cmd_simu_rep));
			if (simuRep < 0) {
				System.err
						.println("simulation replication should be greater than zero");
				System.exit(0);
			}
		}

//simulation polygenic
		if (cl.hasOption(cmd_simu_poly_qt)) {
			simupolyQTFlag = true;
		}

		if (cl.hasOption(cmd_simu_poly_cc)) {
			simupolyCCFlag = true;
		}

		if (cl.hasOption(cmd_poly_loci)) {
			polyLoci = Integer.parseInt(cl.getOptionValue(cmd_poly_loci));
		}

		if(cl.hasOption(cmd_poly_LD)) {
			polyLD = Double.parseDouble(cl.getOptionValue(cmd_poly_LD));
		}

		if(cl.hasOption(cmd_poly_U)) {
			polyU = true;
		}

		if(cl.hasOption(cmd_poly_sample)) {
			polySample = Integer.parseInt(cl.getOptionValue(cmd_poly_sample)); 
		}

		if(cl.hasOption(cmd_poly_case)) {
			polyCS = Integer.parseInt(cl.getOptionValue(cmd_poly_case));
		}
		
		if(cl.hasOption(cmd_poly_k)) {
			polyK = Double.parseDouble(cl.getOptionValue(cmd_poly_k));
		}

		if(cl.hasOption(cmd_poly_hsq)) {
			polyHsq = Double.parseDouble(cl.getOptionValue(cmd_poly_hsq));
		}

//haseman-elston regression
		if(cl.hasOption(cmd_he)) {
			heFlag = true;
		}
		
		if (cl.hasOption(cmd_grm)) {
			StringBuilder sb1 = new StringBuilder(cl.getOptionValue(cmd_grm));
			grm = sb1.append(".grm.gz").toString();
			StringBuilder sb2 = new StringBuilder(cl.getOptionValue(cmd_grm));
			grm_id = sb2.append(".grm.id").toString();
		}

		if (cl.hasOption(cmd_pheno)) {
			pheno = cl.getOptionValue(cmd_pheno);
		}

		if (cl.hasOption(cmd_mpheno)) {
			String[] s = cl.getOptionValue(cmd_mpheno).split(sep);
			mpheno = new int[s.length];
			for (int i = 0; i < s.length; i++) {
				mpheno[i] = Integer.parseInt(s[i]);
			}
		}

		if (cl.hasOption(cmd_reverse)) {
			reverse = true;
		}

		if (cl.hasOption(cmd_he_sd)) {
			heType[he_sd] = true;
			heType[he_ss] = false;
			heType[he_cp] = false;
		}

		if (cl.hasOption(cmd_he_ss)) {
			heType[he_sd] = false;
			heType[he_ss] = true;
			heType[he_cp] = false;
		}

		if (cl.hasOption(cmd_he_cp)) {
			heType[he_sd] = false;
			heType[he_ss] = false;
			heType[he_cp] = true;
		}

		if (cl.hasOption(cmd_k)) {
			k_button = true;
			k = Double.parseDouble(cl.getOptionValue(cmd_k));
		}

		if (cl.hasOption(cmd_na)) {
			na = cl.getOptionValue(cmd_na).split(sep);
		}

		if (cl.hasOption(cmd_out)) {
			out = cl.getOptionValue(cmd_out);
		}

		if (cl.hasOption(cmd_cal_k)) {
			calOption = true;
			cal_k = Double.parseDouble(cl.getOptionValue(cmd_cal_k));
		}

		if (cl.hasOption(cmd_cal_ho)) {
			cal_ho = Double.parseDouble(cl.getOptionValue(cmd_cal_ho));
		}

		if (cl.hasOption(cmd_cal_hl)) {
			cal_hl = Double.parseDouble(cl.getOptionValue(cmd_cal_hl));
		}

		if (cl.hasOption(cmd_cal_cc)) {
			String[] s = cl.getOptionValue(cmd_cal_cc).split(sep);
			cal_cc[0] = Double.parseDouble(s[0]);
			cal_cc[1] = Double.parseDouble(s[1]);
		}

		if (cl.hasOption(cmd_cal_h2_se)) {
			cal_h2_se = Double.parseDouble(cl.getOptionValue(cmd_cal_h2_se));
		}
	}

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}

	public static boolean isNA(String n) {
		boolean f = false;
		for (int i = 0; i < na.length; i++) {
			if(n.compareTo(na[i]) == 0) {
				f = true;
				break;
			}
		}
		return f;
	}
}
