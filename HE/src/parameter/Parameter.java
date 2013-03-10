package parameter;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
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

// singleton implemented in enum way
public enum Parameter {
	INSTANCE;

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

	private final String cmd_realcheck_marker_number = "realcheck_marker_number";
	private final String cmd_realcheck_marker_number_long = "realcheck-marker-number";
	public static int realcheckMarkerNumber = 100;

	private final String cmd_realcheck_snps = "realcheck_snps";
	private final String cmd_realcheck_snps_long = "realcheck-snps";
	public static String realcheckSNPs = null;
	
	private final String cmd_realcheck = "realcheck";
	public static boolean realcheckFlag = false;

	private final String cmd_bfile2 = "bfile2";
	public static String bfile2 = null;
	public static String bedfile2 = null;
	public static String bimfile2 = null;
	public static String famfile2 = null;
	public static boolean bfileOption2 = false;

///////////////merge
	private final String cmd_merge = "merge";
	public static boolean mergeFlag = false;
	private final String cmd_merge_maf_cutoff = "merge_maf_cutoff";
	private final String cmd_merge_maf_cutoff_long = "merge-maf-cutoff";
	public static double merge_maf_cutoff = 0.4;

	private final String cmd_merge_p_cutoff = "merge_p_cutoff";
	private final String cmd_merge_p_cutoff_long = "merge-p-cutoff";
	public static double merge_p_cutoff = 0.05;

	private final String cmd_keep_atgc = "keep_atgc";
	private final String cmd_keep_atgc_long = "keep-atgc";
	public static boolean keepATGCFlag = false;
	
	private final String cmd_remove_Flip = "remove_flip";
	private final String cmd_remove_Flip_long = "remove-flip";
	public static boolean removeFlipFlag = false;

//strand
	private final String cmd_strand = "strand";
	public static boolean strandFlag = false;
	public static String strand_file = null;

//make-predictor
	private final String cmd_make_predictor = "build_predictor";
	private final String cmd_make_predictor_long = "build-predictor";
	public static boolean makePredictorFlag = false;

	private final String cmd_make_predictor2 = "build_predictor2";
	private final String cmd_make_predictor2_long = "build-predictor2";
	public static boolean makePredictor2Flag = false;

	
	private final String cmd_predictor_idx = "predictor_idx";
	private final String cmd_predictor_idx_long = "predictor-idx";
	public static int predictor_idx = 0;

	private final String cmd_predictor_file = "predictor_file";
	private final String cmd_predictor_file_long = "predictor-file";
	public static String predictor_file = null;

	private final String cmd_linear = "linear";
	private final String cmd_logit = "logit";
	public static int LINEAR = 0;
	public static int LOGIT = 1;
	public static int tranFunction = LINEAR;

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

	private final String cmd_simu_realdata = "simu_real_data";
	private final String cmd_simu_realdata_long = "simu-real-data";
	public static boolean simuRealData = false;

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
	public static boolean simupolyQTFlag = false;

	private final String cmd_simu_order = "simu_order";
	private final String cmd_simu_order_long = "simu-order";
	public static boolean simuOrderFlag = false;

	private final String cmd_simu_cc = "simu_cc";
	private final String cmd_simu_cc_long = "simu-cc";
	public static int[] simuCC = { 0, 0 };
	public static boolean simupolyCCFlag = false;

	private final String cmd_simu_k = "simu_k";
	private final String cmd_simu_k_long = "simu-k";
	public static double simuK = 0.1;

	public final static int sm_qt = 0;
	public final static int sm_cc = 1;

	public static boolean[] simuType = { true, false };  //first one for case-control, second one for quantitative

/////////////////simulation polygenic
	private final String cmd_poly_loci = "poly_loci";
	private final String cmd_poly_loci_long = "poly-loci";
	
	public static int polyLoci = 1000;	

	private final String cmd_poly_loci_null = "poly_loci_null";
	private final String cmd_poly_loci_null_long = "poly-loci-null";

	public static int polyLociNull = 0;
	public static int poly_sample_QT = 1000;

	private final String cmd_poly_LD = "poly_ld";
	private final String cmd_poly_LD_long = "poly-ld";
	public static double polyLD = 0;

	private final String cmd_poly_U = "poly_U";
	private final String cmd_poly_U_long = "poly-U";
	public static boolean polyU = false;

	private final String cmd_poly_freq = "poly_freq";
	private final String cmd_poly_freq_long = "poly-freq";
	public static double polyFreq = 0.5;
	
	private final String cmd_poly_effect = "poly_effect";
	private final String cmd_poly_effect_long = "poly-effect";
	public static boolean polyEffectFlag = false;
	public static String polyEffectFile = null;

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
	public static boolean freqFlag = false;
	private final String cmd_geno_freq = "geno_freq";
	private final String cmd_geno_freq_long = "geno-freq";
	public static boolean genoFreqFlag = false;

	public final String cmd_sum_stat_help = "sum_stat_help";
	public final String cmd_sum_stat_help_long = "sum-stat-help";
	public final static int freq = 0;
	public final static int geno_freq = 1;
	public static boolean sumStatFlag = false;
	
	//fst
	private final String cmd_fst = "fst";
	public static boolean fstFlag = false;
	public static String fst_file = null;

//profile
	
	private final String cmd_score = "score";
	private final String cmd_MaCH_Infor = "mach_infor";
	private final String cmd_MaCH_Infor_long = "mach-infor";
	private final String cmd_MaCH_Dosage = "mach_dosage";
	private final String cmd_MaCH_Dosage_long = "mach-dosage";

	private final String cmd_MaCH_Infor_Batch = "mach_infor_batch";
	private final String cmd_MaCH_Infor_Batch_long = "mach-infor-batch";
	private final String cmd_MaCH_Dosage_Batch = "mach_dosage_batch";
	private final String cmd_MaCH_Dosage_Batch_long = "mach-dosage-batch";

	private final String cmd_q_score_file = "q_score_file";
	private final String cmd_q_score_file_long = "q-score-file";
	private final String cmd_q_score_range ="q_score_range";
	private final String cmd_q_score_range_long = "q-score-range";
	
	public static boolean scoreFlag = false;
	public static String scoreFile = null;
	public static String MaCH_Infor = null;
	public static String MaCH_Dosage = null;
	public static String MaCH_Infor_Batch = null;
	public static String MaCH_Dosage_Batch = null;
	public static String q_score_file = null;
	public static String q_score_range_file = null;

//he regression
	private final String cmd_he_sd = "he_sd"; //(y1-y2)^2
	private final String cmd_he_sd_long = "he-sd";
	private final String cmd_he_ss = "he_ss"; //(y1+y2)^2
	private final String cmd_he_ss_long = "he-ss";
	private final String cmd_he_cp = "he_cp"; //y1*y2
	private final String cmd_he_cp_long = "he-cp";
	public static int he_sd = 0;
	public static int he_ss = 1;
	public static int he_cp = 2;
	public static boolean[] heType = {true, false, false};
	public static boolean heFlag = false;
	
	private final String cmd_grm_bin = "grm_bin";
	private final String cmd_grm_bin_long = "grm-bin";
	public static String grm_bin = null;
	public static boolean grm_bin_flag = false;
	
	private final String cmd_grm = "grm";
	public static String grm = null;
	public static String grm_id = null;

	private final String cmd_pheno = "pheno";
	public static String pheno = null;

	private final String cmd_mpheno = "mpheno";
	public static int[] mpheno = {1};
	
	//quantitative covariates
	private final String cmd_qcovar = "qcovar";
	public static String qcovar_file = null;
	private final String cmd_qcovar_num = "qcovar_num";
	private final String cmd_qcovar_num_long = "qcovar-num";
	public static int[] qcovar_num = null;

	//categorical covariates
	private final String cmd_covar = "covar";
	public static String covar_file = null;
	private final String cmd_covar_num = "covar_num";
	private final String cmd_covar_num_long = "covar-num";
	public static int[] covar_num = null;

	private final String cmd_reverse = "reverse";
	public static boolean reverse = false;

	private final String cmd_k = "k";
	public static boolean k_button = false;
	public static double k = 0.01;

	private final String cmd_scale = "scale";
	public static boolean scale = false;

	public static double eh2=1;
	private final String cmd_eh2 = "eh2";
	public static boolean eh2Flag = false;

	private final String cmd_out = "out";
	public static String out = "he";
	
	private final String cmd_perm = "perm";
	public static int perm = 100;
	public static boolean permFlag = false;

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
	public static String[] na = {"-9", "NA", "na", "-Inf", "Inf"};

///////////////////level 1 snp selection
	private final String cmd_chr = "chr";
	public static String[] inchr = null;
	public static String[] exchr = null;
	public static boolean inchrFlag = false;
	public static boolean exchrFlag = false;

	private final String cmd_snps = "snps";
	public static String snpList = null;

//////////////////level 1 individual selection
	private final String cmd_keep = "keep";
	public static String keepFile = null;
	public static boolean keepFlag = false;

///////////////// individual selection start
	private final String cmd_remove = "remove";
	public static String removeFile = null;
	public static boolean removeFlag = false;

///////////////// reference-allele
	private final String cmd_reference_allele = "reference_allele";
	private final String cmd_reference_allele_long = "refernce-allele";
	public static String reference_allele = null;
	
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

	public static String version = "\n"
			+ "******************************************************************\n"
			+ "| HE Jul/18/2011                                                 |\n"
			+ "| (C) 2011 Guo-Bo Chen                                           |\n"
			+ "| v 0.7.7                                                        |\n"			
			+ "| GNU General Public License, v2                                 |\n"
			+ "******************************************************************\n";
	private Options ops = new Options();

	private CommandLineParser parser = new PosixParser();

	private Parameter() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	@SuppressWarnings("static-access")
	public void commandInitial() {

		ops.addOption(OptionBuilder.withDescription("file ").hasArg().create(cmd_file));

		ops.addOption(OptionBuilder.withDescription("bfile ").hasArg().create(cmd_bfile));

//realcheck
		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_upper_long).withDescription("realcheck marker threshold upper bounder").hasArg().create(cmd_realcheck_threshold_upper));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_lower_long).withDescription("realcheck marker threshold lower bounder").hasArg().create(cmd_realcheck_threshold_lower));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_threshold_long).withDescription("realcheck marker threshold").hasArg().create(cmd_realcheck_threshold));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_marker_number_long).withDescription("realcheck marker number").hasArg().create(cmd_realcheck_marker_number));

		ops.addOption(OptionBuilder.withLongOpt(cmd_realcheck_snps_long).withDescription("realcheck snp number").hasArg().create(cmd_realcheck_snps));

		ops.addOption(OptionBuilder.withDescription("realcheck ").create(cmd_realcheck));

		ops.addOption(OptionBuilder.withDescription("bfile2 ").hasArg().create(cmd_bfile2));


//simulation merge
		ops.addOption(OptionBuilder.withDescription("merge ").create(cmd_merge));

		ops.addOption(OptionBuilder.withLongOpt(cmd_merge_maf_cutoff_long).withDescription("merge maf cutoff").hasArg().create(cmd_merge_maf_cutoff));

		ops.addOption(OptionBuilder.withLongOpt(cmd_merge_p_cutoff_long).withDescription("merge p cutoff").hasArg().create(cmd_merge_p_cutoff));

		ops.addOption(OptionBuilder.withLongOpt(cmd_keep_atgc_long).withDescription("remove A/T and G/C loci").create(cmd_keep_atgc));

		ops.addOption(OptionBuilder.withLongOpt(cmd_remove_Flip_long).withDescription("remove flipped loci").create(cmd_remove_Flip));
//make predictor
		ops.addOption(OptionBuilder.withLongOpt(cmd_make_predictor_long).withDescription("make predictor").create(cmd_make_predictor));

		ops.addOption(OptionBuilder.withLongOpt(cmd_make_predictor2_long).withDescription("make predictor2").create(cmd_make_predictor2));

		ops.addOption(OptionBuilder.withLongOpt(cmd_predictor_idx_long).withDescription("predictor index").hasArg().create(cmd_predictor_idx));

		ops.addOption(OptionBuilder.withLongOpt(cmd_predictor_file_long).withDescription("predictor file").hasArg().create(cmd_predictor_file));

		ops.addOption(OptionBuilder.withDescription("linear").create(cmd_linear));

		ops.addOption(OptionBuilder.withDescription("logit").create(cmd_logit));

//simulation nuclear fam
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_long).withDescription("simulation nuclear family ").create(cmd_simu_fam));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_size_long).withDescription("simulation nuclear family size ").hasArg().create(cmd_simu_fam_size));
		
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_fam_marker_long).withDescription("simulation number for nuclear family ").hasArg().create(cmd_simu_fam_marker));		
		
//simulation real data
		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_realdata_long).withDescription("gwas simulations ").hasArg().create(cmd_simu_realdata));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_seed_long).withDescription("gwas simulation seed ").hasArg().create(cmd_simu_seed));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_rep_long).withDescription("gwas simulation replication ").hasArg().create(cmd_simu_rep));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_casual_loci_long).withDescription("gwas simulation casual loci ").hasArg().create(cmd_simu_casual_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_rnd_casual_loci_long).withDescription("gwas simulation casual loci number ").hasArg().create(cmd_simu_rnd_casual_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_hsq_long).withDescription("gwas simulation heritability ").hasArg().create(cmd_simu_hsq));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_qt_long).withDescription("gwas simulate quantitative traits ").hasArg().create(cmd_simu_qt));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_order_long).withDescription("order SNP effects ascendingly ").create(cmd_simu_order));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_cc_long).withDescription("gwas simulate case-control ").hasArgs(2).create(cmd_simu_cc));

		ops.addOption(OptionBuilder.withLongOpt(cmd_simu_k_long).withDescription("gwas prevalence of the binary trait ").hasArg().create(cmd_simu_k));

//nontransmitted
		ops.addOption(OptionBuilder.withDescription("nontransmitted ").create(cmd_nontrans));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_cases_long).withDescription("nontransmitted filter cases ").create(cmd_nontrans_cases));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_controls_long).withDescription("nontransmitted filter controls ").create(cmd_nontrans_controls));

		ops.addOption(OptionBuilder.withLongOpt(cmd_nontrans_seed_long).withDescription("gwas prevalence of the binary trait ").hasArg().create(cmd_nontrans_seed));
		
//simulation polygenic model

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_loci_long).withDescription("number of polygenic loci, defualt= " + polyLoci).hasArg().create(cmd_poly_loci));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_loci_null_long).withDescription("number of null polygenic loci, defualt= " + polyLoci).hasArg().create(cmd_poly_loci_null));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_LD_long).withDescription("LD (correlation), defualt= " + polyLD).hasArg().create(cmd_poly_LD));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_U_long).withDescription("polygenic model has Uniform Effect? " + polyU).create(cmd_poly_U));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_freq_long).withDescription("minor allele frequency for polygenic model? " + polyFreq).hasArg().create(cmd_poly_freq));

		ops.addOption(OptionBuilder.withLongOpt(cmd_poly_effect_long).withDescription("effect for polygenic model? " + polyEffectFile).hasArg().create(cmd_poly_effect));

//pop stat
		ops.addOption(OptionBuilder.withDescription("calculate MAF frequency ").create(cmd_freq));

		ops.addOption(OptionBuilder.withLongOpt(cmd_geno_freq_long).withDescription("calculate genotype frequency ").create(cmd_geno_freq));

		ops.addOption(OptionBuilder.withDescription("calculate fst ").hasArg().create(cmd_fst));

//snp selection
		ops.addOption(OptionBuilder.withDescription("select chromosomes").hasArgs().create(cmd_chr));
		ops.addOption(OptionBuilder.withDescription("select snps").hasArgs().create(cmd_snps));

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

		ops.addOption(OptionBuilder.withLongOpt(cmd_reference_allele_long).withDescription("set reference allele ").hasArg().create(cmd_reference_allele));

		ops.addOption(OptionBuilder.withLongOpt(cmd_make_bed_long).withDescription("make bed ").create(cmd_make_bed));

		ops.addOption(OptionBuilder.withDescription("solve strand ").hasArg().create(cmd_strand));

//profile
		ops.addOption(OptionBuilder.withDescription("profile score").hasArg().create(cmd_score));
		
		ops.addOption(OptionBuilder.withLongOpt(cmd_MaCH_Dosage_long).withDescription("Mach dosage file").hasArg().create(cmd_MaCH_Dosage));

		ops.addOption(OptionBuilder.withLongOpt(cmd_MaCH_Infor_long).withDescription("Mach dosage information file").hasArg().create(cmd_MaCH_Infor));

		ops.addOption(OptionBuilder.withLongOpt(cmd_MaCH_Dosage_Batch_long).withDescription("Mach dosage batch file").hasArg().create(cmd_MaCH_Dosage_Batch));

		ops.addOption(OptionBuilder.withLongOpt(cmd_MaCH_Infor_Batch_long).withDescription("Mach dosage information batch file").hasArg().create(cmd_MaCH_Infor_Batch));

		ops.addOption(OptionBuilder.withLongOpt(cmd_q_score_file_long).withDescription("q score file").hasArg().create(cmd_q_score_file));

		ops.addOption(OptionBuilder.withLongOpt(cmd_q_score_range_long).withDescription("q score range").hasArg().create(cmd_q_score_range));

//haseman-elston regression
		ops.addOption(OptionBuilder.withDescription("h2 ").hasArg().create(cmd_eh2));

		ops.addOption(OptionBuilder.withLongOpt(cmd_grm_bin_long).withDescription("grm binary format").hasArg().create(cmd_grm_bin));
		
		ops.addOption(OptionBuilder.withDescription("grm ").hasArg().create(cmd_grm));

		ops.addOption(OptionBuilder.withDescription("pheno ").hasArg().create(cmd_pheno));

		ops.addOption(OptionBuilder.withLongOpt(cmd_mpheno).withDescription("pheno number " + cmd_mpheno).hasArg().withArgName("index").create(cmd_mpheno));

		ops.addOption(OptionBuilder.withDescription("covariate file").hasArg().create(cmd_covar));

		ops.addOption(OptionBuilder.withLongOpt(cmd_covar_num_long).withDescription("covariate index").hasArg().create(cmd_covar_num));

		ops.addOption(OptionBuilder.withDescription("quantitative covariate file").hasArg().create(cmd_qcovar));

		ops.addOption(OptionBuilder.withLongOpt(cmd_qcovar_num_long).withDescription("quantitative covariate index").hasArg().create(cmd_qcovar_num));

		ops.addOption(OptionBuilder.withDescription("reverse ").create(cmd_reverse));

		ops.addOption(OptionBuilder.withDescription("standardise the phenotype").create(cmd_scale));

		ops.addOption(OptionBuilder.withDescription("perm ").hasArg().create(cmd_perm));

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

	public void commandListener(String[] args) {
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

		if (cl.hasOption(cmd_realcheck_marker_number)) {
			realcheckMarkerNumber = Integer.parseInt(cl.getOptionValue(cmd_realcheck_marker_number));
			if (realcheckMarkerNumber < 0) {
				System.err.println("realcheck marker number should be greater than 0");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_realcheck_snps)) {
			realcheckSNPs = cl.getOptionValue(cmd_realcheck_snps);
			exists(realcheckSNPs);
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
		
		if (cl.hasOption(cmd_snps)) {
			snpList = cl.getOptionValue(cmd_snps);
			exists(snpList);
		}
		
//individual selection 1 keep
		if (cl.hasOption(cmd_keep)) {
			keepFile = cl.getOptionValue(cmd_keep);
			exists(keepFile);
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
			exists(removeFile);
			removeFlag = true;
		}

/////// set reference
		if (cl.hasOption(cmd_reference_allele)) {
			reference_allele = cl.getOptionValue(cmd_reference_allele);
			exists(reference_allele);
		}
		
		if (cl.hasOption(cmd_strand)) {
			strand_file = cl.getOptionValue(cmd_strand);
			exists(strand_file);
			strandFlag = true;
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
			freqFlag = true;
		}
		if (cl.hasOption(cmd_geno_freq)) {
			sumStatFlag = true;
			genoFreqFlag = true;
		}
		if (cl.hasOption(cmd_fst)) {
			sumStatFlag = true;
			fstFlag = true;
			fst_file = cl.getOptionValue(cmd_fst);
			exists(fst_file);
		}

//merge 
		if (cl.hasOption(cmd_merge)) {
			mergeFlag = true;
		}
		
		if (cl.hasOption(cmd_merge_maf_cutoff)) {
			merge_maf_cutoff = Double.parseDouble(cl.getOptionValue(cmd_merge_maf_cutoff));
			if (merge_maf_cutoff > 0.5 || merge_maf_cutoff < 0) {
				System.err.println("merger maf cutoff should be between 0 and 0.5");
				System.exit(0);
			}
		}		

		if (cl.hasOption(cmd_merge_p_cutoff)) {
			merge_p_cutoff = Double.parseDouble(cl.getOptionValue(cmd_merge_p_cutoff));
			if (merge_p_cutoff < 0) {
				System.err.println("merger p cutoff should be between 0 and 1.");
				System.exit(0);
			}
		}		

		if (cl.hasOption(cmd_keep_atgc)) {
			keepATGCFlag = true;
		}

		if (cl.hasOption(cmd_remove_Flip)) {
			removeFlipFlag = true;
		}

//make predictor
		if (cl.hasOption(cmd_make_predictor)) {
			makePredictorFlag = true;
		}

		if (cl.hasOption(cmd_make_predictor2)) {
			makePredictor2Flag = true;
		}

		if (cl.hasOption(cmd_predictor_idx)) {
			predictor_idx = Integer.parseInt(cl.getOptionValue(cmd_predictor_idx));
			predictor_idx--;
			if (predictor_idx < 0) {
				System.err.println("predictor-idx should bigger than 0");
				System.exit(0);
			}
		}

		if (cl.hasOption(cmd_predictor_file)) {
			predictor_file = cl.getOptionValue(cmd_predictor_file);
			exists(predictor_file);
		}

		if (cl.hasOption(cmd_linear)) {
			tranFunction = LINEAR;
		}

		if (cl.hasOption(cmd_logit)) {
			tranFunction = LOGIT;
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

		if (cl.hasOption(cmd_simu_realdata)) {
			simuRealData = true;
		}
		
		if (cl.hasOption(cmd_simu_qt)) {
			simuType[sm_qt] = false;
			simuType[sm_cc] = true;
			simupolyQTFlag = true;
			
			poly_sample_QT = Integer.parseInt(cl.getOptionValue(cmd_simu_qt));
		}

		if (cl.hasOption(cmd_simu_cc)) {
			simuType[sm_qt] = true;
			simuType[sm_cc] = false;
			
			String[] s = cl.getOptionValues(cmd_simu_cc);
			simuCC[0] = Integer.parseInt(s[0]);
			simuCC[1] = Integer.parseInt(s[1]);
			simupolyCCFlag = true;
		}

		if (cl.hasOption(cmd_simu_order)) {
			simuOrderFlag = true;
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

		if (cl.hasOption(cmd_poly_loci)) {
			polyLoci = Integer.parseInt(cl.getOptionValue(cmd_poly_loci));
		}

		if (cl.hasOption(cmd_poly_loci_null)) {
			polyLociNull = Integer.parseInt(cl.getOptionValue(cmd_poly_loci_null));
		}

		if (cl.hasOption(cmd_poly_LD)) {
			polyLD = Double.parseDouble(cl.getOptionValue(cmd_poly_LD));
		}

		if (cl.hasOption(cmd_poly_U)) {
			polyU = true;
		}

		if (cl.hasOption(cmd_poly_freq)) {
			polyFreq = Double.parseDouble(cl.getOptionValue(cmd_poly_freq));
		}

		if (cl.hasOption(cmd_poly_effect)) {
			polyEffectFlag = true;
			polyEffectFile = cl.getOptionValue(cmd_poly_effect);
			exists(polyEffectFile);
		}

//profile
		if (cl.hasOption(cmd_score)) {
			scoreFlag = true;
			scoreFile = cl.getOptionValue(cmd_score);
			exists(scoreFile);
		}

		if (cl.hasOption(cmd_MaCH_Infor)) {
			MaCH_Infor = cl.getOptionValue(cmd_MaCH_Infor);
			exists(MaCH_Infor);
		}
		
		if (cl.hasOption(cmd_MaCH_Dosage)) {
			MaCH_Dosage = cl.getOptionValue(cmd_MaCH_Dosage);
			exists(MaCH_Dosage);
		}

		if (cl.hasOption(cmd_MaCH_Infor_Batch)) {
			MaCH_Infor_Batch = cl.getOptionValue(cmd_MaCH_Infor_Batch);
			exists(MaCH_Infor_Batch);
		}

		if (cl.hasOption(cmd_MaCH_Dosage_Batch)) {
			MaCH_Dosage_Batch = cl.getOptionValue(cmd_MaCH_Dosage_Batch);
			exists(MaCH_Dosage_Batch);
		}

		if (cl.hasOption(cmd_q_score_file)) {
			q_score_file = cl.getOptionValue(cmd_q_score_file);
			exists(q_score_file);
		}

		if (cl.hasOption(cmd_q_score_range)) {
			q_score_range_file = cl.getOptionValue(cmd_q_score_range);
			exists(q_score_range_file);
		}

//haseman-elston regression

		if (cl.hasOption(cmd_eh2)) {
			eh2Flag = true;
			eh2 = Double.parseDouble(cl.getOptionValue(cmd_eh2));
		}

		if (cl.hasOption(cmd_grm)) {
			StringBuilder sb1 = new StringBuilder(cl.getOptionValue(cmd_grm));
			grm = sb1.append(".grm.gz").toString();
			StringBuilder sb2 = new StringBuilder(cl.getOptionValue(cmd_grm));
			grm_id = sb2.append(".grm.id").toString();
			grm_bin_flag = false;
		}

		if (cl.hasOption(cmd_grm_bin)) {
			StringBuilder sb1 = new StringBuilder(cl.getOptionValue(cmd_grm_bin));
			grm_bin = sb1.append(".grm.bin").toString();
			StringBuilder sb2 = new StringBuilder(cl.getOptionValue(cmd_grm_bin));
			grm_id = sb2.append(".grm.id").toString();
			grm_bin_flag = true;
		}

		if (cl.hasOption(cmd_pheno)) {
			pheno = cl.getOptionValue(cmd_pheno);
		}

		if (cl.hasOption(cmd_mpheno)) {
			String[] s = cl.getOptionValue(cmd_mpheno).split(",");
			mpheno = new int[s.length];
			for (int i = 0; i < s.length; i++) {
				mpheno[i] = Integer.parseInt(s[i]);
			}
		}

		if (cl.hasOption(cmd_covar)) {
			covar_file = cl.getOptionValue(cmd_covar);
			exists(covar_file);
		}

		if (cl.hasOption(cmd_covar_num)) {
			String[] p = cl.getOptionValue(cmd_covar_num).split(",");
			HashSet<Integer> idx = NewIt.newHashSet();
			for (int i = 0, len = p.length; i < len; i++) {
				if (p[i].contains("-")) {
					String[] pp = p[i].split("-");
					if (pp.length != 2) {
						System.err.println("bad parameter for option --" + cmd_covar_num_long + ": " + p[i] +".");
						Test.LOG.append("bad parameter for option --" + cmd_covar_num_long + ": " + p[i] +".\n");
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
			covar_num = new int[idx.size()];
			int c = 0;
			for (Iterator<Integer> e = idx.iterator(); e.hasNext();) {
				covar_num[c] = e.next().intValue();
				if (covar_num[c] < 0) {
					System.err.println("bad parameter for option --" + cmd_covar_num_long + ": " + covar_num[c] +".");
					Test.LOG.append("bad parameter for option --" + cmd_covar_num_long + ": " + covar_num[c] +".\n");
					Test.printLog();
					System.exit(0);
				}
				c++;
			}
		}

		if (cl.hasOption(cmd_qcovar)) {
			qcovar_file = cl.getOptionValue(cmd_qcovar);
			exists(qcovar_file);
		}

		if (cl.hasOption(cmd_qcovar_num)) {
			String[] p = cl.getOptionValue(cmd_qcovar_num).split(",");
			HashSet<Integer> idx = NewIt.newHashSet();
			for (int i = 0, len = p.length; i < len; i++) {
				if (p[i].contains("-")) {
					String[] pp = p[i].split("-");
					if (pp.length != 2) {
						System.err.println("bad parameter for option --" + cmd_qcovar_num_long + ": " + p[i] +".");
						Test.LOG.append("bad parameter for option --" + cmd_qcovar_num_long + ": " + p[i] +".\n");
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
			qcovar_num = new int[idx.size()];
			int c = 0;
			for (Iterator<Integer> e = idx.iterator(); e.hasNext();) {
				qcovar_num[c] = e.next().intValue();
				if (qcovar_num[c] < 0) {
					System.err.println("bad parameter for option --" + cmd_qcovar_num_long + ": " + qcovar_num[c] +".");
					Test.LOG.append("bad parameter for option --" + cmd_qcovar_num_long + ": " + qcovar_num[c] +".\n");
					Test.printLog();
					System.exit(0);
				}
				c++;
			}
		}

		if (cl.hasOption(cmd_reverse)) {
			reverse = true;
		}

		if (cl.hasOption(cmd_scale)) {
			scale = true;
		}

		if (cl.hasOption(cmd_perm)) {
			permFlag = true;
			perm = Integer.parseInt(cl.getOptionValue(cmd_perm));
		}

		if (cl.hasOption(cmd_he_sd)) {
			heFlag = true;
			heType[he_sd] = true;
			heType[he_ss] = false;
			heType[he_cp] = false;
		}

		if (cl.hasOption(cmd_he_ss)) {
			heFlag = true;
			heType[he_sd] = false;
			heType[he_ss] = true;
			heType[he_cp] = false;
		}

		if (cl.hasOption(cmd_he_cp)) {
			heFlag = true;
			heType[he_sd] = false;
			heType[he_ss] = false;
			heType[he_cp] = true;
		}

		if (cl.hasOption(cmd_k)) {
			k_button = true;
			k = Double.parseDouble(cl.getOptionValue(cmd_k));
		}

		if (cl.hasOption(cmd_na)) {
			na = cl.getOptionValue(cmd_na).split(",");
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
			String[] s = cl.getOptionValue(cmd_cal_cc).split(",");
			cal_cc[0] = Double.parseDouble(s[0]);
			cal_cc[1] = Double.parseDouble(s[1]);
		}

		if (cl.hasOption(cmd_cal_h2_se)) {
			cal_h2_se = Double.parseDouble(cl.getOptionValue(cmd_cal_h2_se));
		}
	}

	public static void main(String[] args) {
		Parameter.INSTANCE.commandListener(args);
		System.out.println(Parameter.INSTANCE);
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
	
	private void exists(String file) {
		File f = new File(file);
		if (!f.exists()) {
			System.err.println("could not open " + file + ".");
			System.exit(0);
		}
	}
}
