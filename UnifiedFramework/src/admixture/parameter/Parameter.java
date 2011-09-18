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

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Parameter {

	private final String incommand_separator = ",";
	private final String cmd_missing_allele = "missingallele";
	public static String missing_allele = "0";

	private final String cmd_missing_phenotype = "missingphenotype";
	public static String missing_phenotype = "-9";

	private final String cmd_missing_genotype = "missinggenotype";
	public static String missing_genotype = "0";

	private final String cmd_status_shift = "1";
	public static int status_shift = 0;

	private final String cmd_mode = "md";// u, f (sii), pi(ajhg2008);
	public static String mode = "u";

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
	private final String cmd_tfile = "tfile";
	public static boolean tfileFlag = false;
	private final String cmd_tped = "tped";
	public static String tped = null;
	private final String cmd_tfam = "tfam";
	public static String tfam = null;
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
	
	private final String cmd_method = "model";
	public int linkfunction = 0;
	private String[] lf = new String[] { "0 for linear regression", "1 for logistic regression" };
	// phenotype set end

	// individual selection start
	private final String cmd_ex_fam = "exfam";
	public static String[] ex_family = null;
	private final String cmd_ex_fam_file = "exfamfile";
	public static boolean exfamFlag = false;

	private final String cmd_ex_ind = "exind";
	public static String[][] ex_ind = null;
	private final String cmd_ex_ind_file = "exindfile";
	public static boolean exindFlag = false;

	private final String cmd_filter_male = "filtermale";
	public static boolean filter_maleFlag = false;

	private final String cmd_filter_female = "filterfemale";
	public static boolean filter_femaleFlag = false;

	private final String cmd_ex_nosex = "exnosex";
	public static boolean ex_nosexFlag = false;
	// end individual filter

	// snp selection
	private final String cmd_chr = "chr";
	public static String[] chr = null;
	public static boolean chrFlag = false;

	private final String cmd_snpwindow = "snpwindow";
	public static String[] snpwindow = null;
	public static double[][] snp_window = null;
	public static boolean snpwindowFlag = false;

	private final String cmd_snp = "snp";
	private final String cmd_snp_f = "snpfile";
	public static boolean snpFlag = false;
	public static String[] includesnp = null;
	public static String[] excludesnp = null;

	public static String[] insnpPair = null;
	public static String[] exsnpPair = null;
	public static boolean snpPairFlag = false;

	   //this set only used when the option x is specified;
	public static String[] xincludesnp = null;
	public static String[] xinsnpPair = null;

	public static boolean xsnpFlag = false;
	public static boolean xsnpPairFlag = false;
	   //end it

	private final String cmd_bgsnp = "bgsnp";
	public static boolean bgsnpFlag = false;
	public static String[] bgsnp = null;

	private final String cmd_x = "x";
	public static boolean x = false;
	// snp selection end

	//soft snp selection 
	private final String cmd_maf = "maf";
	public static double maf = -1;
	public static boolean mafFlag = false;
	
	private final String cmd_geno = "geno";
	public static double geno = 2;
	public static boolean genoFlag = false;
	
	private final String cmd_header = "header";
	public static boolean header = false;

	//sampling & partitioning start
	private final String cmd_thin = "thin";
	public static double thin = 1.0;
	
	private final String cmd_slice = "slice";
	public static int sliceN = 1;
	public static int slice = 1;
	//sampling & partitioning end
	
	//mdr options start
	private final String cmd_cv = "cv";
	public static int cv = 5;
	
	private final String cmd_order = "order";
	public static int order = 1;

	private final String cmd_seed = "seed";
	public static int seed = 2011;
	//mdr option end
	
	
	private final String cmd_perm = "permut";
	public int permutation = 100;

	private final String cmd_perm_scheme = "ps";
	public boolean permu_scheme = false;

	private final String cmd_perm_fam = "pf";
	public boolean permu_fam = false;

	private final String cmd_unrelated_only = "ur";
	public boolean unrelated_only = false;

	private final String cmd_simu = "simu";
	public int simu = 1;

	private final String cmd_help = "help";
	public boolean help = false;

	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();

	public Parameter() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	public void commandInitial() {
		ops.addOption(OptionBuilder.withDescription("u (default) for the unified framework and f for using sibs only.").hasArg().create(cmd_mode));

		ops.addOption(OptionBuilder.withDescription("specify the .ped and .map files").hasArg().create(cmd_file));
		ops.addOption(OptionBuilder.withDescription("specify the .ped file.").hasArg().create(cmd_ped));
		ops.addOption(OptionBuilder.withDescription("specify the .map file.").hasArg().create(cmd_map));

		ops.addOption(OptionBuilder.withDescription("specify the .bed, .bim and .fam files.").hasArg().create(cmd_bfile));
		ops.addOption(OptionBuilder.withDescription("specify the .bed file.").hasArg().create(cmd_bed));
		ops.addOption(OptionBuilder.withDescription("specify the .bim file.").hasArg().create(cmd_bim));
		ops.addOption(OptionBuilder.withDescription("specify the .fam file.").hasArg().create(cmd_fam));

		ops.addOption(OptionBuilder.withDescription("specify the .tped and .tfam files").hasArg().create(cmd_tfile));
		ops.addOption(OptionBuilder.withDescription("specify the .tped file.").hasArg().create(cmd_tped));
		ops.addOption(OptionBuilder.withDescription("specify the .tfam file.").hasArg().create(cmd_tfam));

		ops.addOption(OptionBuilder.withDescription("specify the phenotype file.").hasArg().create(cmd_pheno));
		ops.addOption(OptionBuilder.withDescription("specify 1 or more covariates by number.").hasArgs().create(cmd_covar));
		ops.addOption(OptionBuilder.withDescription("specify 1 or more covariates by name.").hasArgs().create(cmd_covar_name));

		ops.addOption(OptionBuilder.withDescription("specify the window size for a snp").hasArgs().create(cmd_snpwindow));

		ops.addOption(OptionBuilder.withDescription("include snps in detecting interaction").hasArgs().create(cmd_snp));
		ops.addOption(OptionBuilder.withDescription("specify the background snp").hasArgs().create(cmd_bgsnp));
		ops.addOption(OptionBuilder.withDescription("specify the file containing included snps when detecting interaction").hasArg().create(
				cmd_snp_f));

		ops.addOption(OptionBuilder.withDescription("specify interacting snps").create(cmd_x));

		ops.addOption(OptionBuilder.withDescription("specify excluded families").hasArgs().create(cmd_ex_fam));
		ops.addOption(OptionBuilder.withDescription("specify the file containing excluded family ids").hasArg().create(cmd_ex_fam_file));
		ops.addOption(OptionBuilder.withDescription("specify excluded individuals").hasArgs().create(cmd_ex_ind));
		ops.addOption(OptionBuilder.withDescription("specify the file containing excluded individual ids").hasArg().create(cmd_ex_ind_file));
		ops.addOption(OptionBuilder.withDescription("include males only").create(cmd_filter_male));
		ops.addOption(OptionBuilder.withDescription("include females only").create(cmd_filter_female));
		ops.addOption(OptionBuilder.withDescription("include unknown sex ").create(cmd_ex_nosex));

		ops.addOption(OptionBuilder.withDescription("select chromosomes").hasArgs().create(cmd_chr));
		ops.addOption(OptionBuilder.withDescription("specify response by number.").hasArg().create(cmd_res_number));
		ops.addOption(OptionBuilder.withDescription("specify response by name.").hasArg().create(cmd_res_name));

		ops.addOption(OptionBuilder.withDescription(
				"method for adjustment of the phenotype, 0 (default) for linear regression, 1 for logistic regression.").hasArg().create(cmd_method));
		ops.addOption(OptionBuilder.withDescription("fold of cross-validation, and default is 5.").hasArg().create(cmd_cv));
		ops.addOption(OptionBuilder.withDescription("specify the order of interaction.").hasArg().create(cmd_order));
		ops.addOption(OptionBuilder.withDescription("specify the random sample fraction.").hasArg().create(cmd_thin));
		ops.addOption(OptionBuilder.withDescription("specify partition of the searching space.").hasArg().create(cmd_slice));
		ops.addOption(OptionBuilder.withDescription("specify the minor allele frequency for inclusion.").hasArg().create(cmd_maf));
		ops.addOption(OptionBuilder.withDescription("specify missing genotype rate for inclusion.").hasArg().create(cmd_geno));
		ops.addOption(OptionBuilder.withDescription("seed for the algorithms").hasArg().create(cmd_seed));

		ops.addOption(OptionBuilder.withDescription("replication for permutation.  Default is 100.").hasArg().create(cmd_perm));
		ops.addOption(OptionBuilder.withDescription("only sibs are exchangeable when this option is turned on").create(cmd_perm_scheme));
		ops.addOption(OptionBuilder.withDescription(
				"hierachical permutation for families that founders are exchangeable within family and sibs are exchangeable within family.").create(
				cmd_perm_fam));
		ops.addOption(OptionBuilder.withDescription("use unrelated indivuduals only, if '--md' is specified.").create(cmd_unrelated_only));
		ops.addOption(OptionBuilder.withDescription("replications for simulation, and this parameter is for simulation only").hasArg().create(
				cmd_simu));
		ops.addOption(OptionBuilder.withDescription("missing phenotype, default 99").hasArg().create(cmd_missing_phenotype));
		ops.addOption(OptionBuilder.withDescription("missing genotype, default 00").hasArg().create(cmd_missing_genotype));
		ops.addOption(OptionBuilder.withDescription("missing allele, default 0").hasArg().create(cmd_missing_allele));
		ops.addOption(OptionBuilder.withDescription("use this option if status was coded as 1 (unaffected)/2 (affected).").create(cmd_status_shift));
		ops.addOption(OptionBuilder.withDescription("help manual.").create(cmd_help));
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
		if (cl.hasOption(cmd_mode)) {
			mode = cl.getOptionValue(cmd_mode);
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
		if (cl.hasOption(cmd_tfile)) {
			StringBuffer sb1 = new StringBuffer();
			StringBuffer sb2 = new StringBuffer();
			sb1.append(cl.getOptionValue(cmd_tfile));
			sb1.append(".tped");

			sb2.append(cl.getOptionValue(cmd_tfile));
			sb2.append(".tfam");

			tped = sb1.toString();
			tfam = sb2.toString();
		}
		if (cl.hasOption(cmd_tped)) {
			tped = cl.getOptionValue(cmd_tped);
		}
		if (cl.hasOption(cmd_tfam)) {
			tfam = cl.getOptionValue(cmd_tfam);
		}
		if (tped != null && tfam != null) {
			File fped = new File(tped);
			if (!fped.exists()) {
				throw new IllegalArgumentException("could not open " + tped);
			}
			File ffam = new File(tfam);
			if (!ffam.exists()) {
				throw new IllegalArgumentException("could not open " + tfam);
			}
			tfileFlag = true;
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
						System.err.println("unknow option value " + p[i] + "\n");
						System.exit(1);
					}
					for (int j = Integer.parseInt(pp[0]); j <= Integer.parseInt(pp[1]); j++) {
						idx.add(new Integer(j));
					}
				} else {
					idx.add(new Integer(Integer.parseInt(p[i])));
				}
			}
			predictor = new int[idx.size()];
			int c = 0;
			for (Iterator<Integer> e = idx.iterator(); e.hasNext();)
				predictor[c++] = e.next().intValue();
		}

		if (cl.hasOption(cmd_covar_name)) {
			HashSet<String> cn = NewIt.newHashSet();
			String[] p = cl.getOptionValues(cmd_covar_name);
			for (int i = 0; i < p.length; i++) {
				if (p[i].contains("-")) {
					String[] pp = predictor_name[i].split("-");
					if (pp.length != 2) {
						throw new IllegalArgumentException("unknown parameter " + predictor_name[i]);
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
								throw new IllegalArgumentException("bad parameter " + snp[i]);
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
								throw new IllegalArgumentException("bad parameter " + snp[i]);
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
						throw new IllegalArgumentException("bad parameter " + snp[i]);
					} else {
						if (snp[i].contains("-")) {
							String[] s = snp[i].split("-");
							if (s.length != 2) {
								throw new IllegalArgumentException("bad parameter " + snp[i]);
							}
							insnppair.add(s[0]);
							insnppair.add(s[1]);
						} else {
							insnp.add(snp[i]);
						}
					}
				}
				if (insnp.size() > 0) {
					xincludesnp = (String[]) insnp.toArray(new String[0]);
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
			for(int i = 0; i < bg.length; i++) {
				bgSet.add(bg[i]);
			}
			if(bgSet.size() != bg.length) {
				throw new IllegalArgumentException("bad parameter for bgsnp");
			}
			bgsnp = cl.getOptionValues(cmd_bgsnp);
			bgsnpFlag = true;
		}

		if (cl.hasOption(cmd_snp_f)) {
			if (!x) {
				String snps_file = cl.getOptionValue(cmd_snp_f);

				File f = new File(snps_file);
				if (!f.exists()) {
					throw new IllegalArgumentException("could not find " + snps_file);
				}
				BufferedReader reader = null;
				try {
					reader = new BufferedReader(new FileReader(f));
				} catch (IOException E) {
					throw new IllegalArgumentException("could not open snps file " + snps_file);
				}

				ArrayList<String> snp = NewIt.newArrayList();
				String line = null;
				try {
					while ((line = reader.readLine()) != null) {
						snp.add(line);
					}
					reader.close();
				} catch (IOException E) {
					throw new IllegalArgumentException("bad lines in " + snps_file);
				}
				if (snp.size() > 0) {
					ArrayList<String> insnp = NewIt.newArrayList();
					ArrayList<String> exsnp = NewIt.newArrayList();
					for (int i = 0; i < snp.size(); i++) {
						if (snp.get(i).startsWith("-")) {
							exsnp.add(snp.get(i).substring(1, snp.get(i).length()));
						} else {
							insnp.add(snp.get(i));
						}
					}
					includesnp = (String[]) insnp.toArray(new String[0]);
					excludesnp = (String[]) exsnp.toArray(new String[0]);
				}
			} else {
				
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

		if (cl.hasOption(cmd_ex_ind)) {
			String[] ind = cl.getOptionValues(cmd_ex_ind);
			ex_ind = new String[2][ind.length];
			for (int i = 0; i < ind.length / 2; i++) {
				String[] Ind = ind[i].split(incommand_separator);
				ex_ind[0][i] = Ind[0];
				ex_ind[1][i] = Ind[1];
			}
			exindFlag = true;
		}
		if (cl.hasOption(cmd_ex_ind_file)) {
			String file = cl.getOptionValue(cmd_ex_ind_file);
			File find = new File(file);
			if (!find.exists()) {
				throw new IllegalArgumentException("could not open file: " + file);
			}
			ArrayList<String> fid = NewIt.newArrayList();
			ArrayList<String> iid = NewIt.newArrayList();

			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(new File(file)));
			} catch (IOException E) {
				throw new IllegalArgumentException("failed in reading " + file);
			}

			String line;
			try {
				while ((line = reader.readLine()) != null) {
					String[] l = line.split(incommand_separator);
					if (l.length < 2) {
						continue;
					}
					fid.add(l[0]);
					iid.add(l[1]);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (fid.size() > 0) {
				exindFlag = true;
				ex_ind = new String[2][fid.size()];
				for (int i = 0; i < ex_ind.length; i++) {
					ex_ind[0][i] = fid.get(i);
					ex_ind[1][i] = iid.get(i);
				}
			} else {
				throw new IllegalArgumentException("bad lines in " + file);
			}
		}
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
			chr = cl.getOptionValues(cmd_chr);
			HashSet<String> chrSet = NewIt.newHashSet();
			for (int i = 0; i < chr.length; i++) {
				chrSet.add(chr[i]);
			}
			if (chr.length != chrSet.size()) {
				throw new IllegalArgumentException("bad parameter for --chr");
			}
			chrFlag = true;
		}

		if (cl.hasOption(cmd_snpwindow)) {
			String[] s = cl.getOptionValues("snpwindow");
			snpwindow = new String[s.length];
			snp_window = new double[s.length][2];
			for (int i = 0; i < s.length; i++) {
				String[] ss = s[i].split(incommand_separator);
				if(ss.length != 3) {
					throw new IllegalArgumentException("bad parameter for --snpwindow: " + s[i]);
				}
				snpwindow[i] = ss[0];
				snp_window[i][0] = Double.parseDouble(ss[1]) * -1000;
				snp_window[i][1] = Double.parseDouble(ss[2]) * 1000;
			}
			snpwindowFlag = true;
		}
		
		if (cl.hasOption(cmd_maf)) {
			maf = Double.parseDouble(cl.getOptionValue(cmd_maf));
			if (maf < 0) {
				throw new IllegalArgumentException("bad parameter for --maf: " + maf);
			}
			mafFlag = true;
		}
		
		if (cl.hasOption(cmd_geno)) {
			geno = Double.parseDouble(cl.getOptionValue(cmd_geno));
			if (geno < 0) {
				throw new IllegalArgumentException("bad parameter for --geno: " + geno);
			}
			genoFlag = true;
		}

		if (cl.hasOption(cmd_header)) {
			header = true;
		}
		// if (cl.hasOption(cmd_topN)) {
		// topN = Integer.parseInt(cl.getOptionValue(cmd_topN));
		// }
		if (cl.hasOption(cmd_res_number)) {
			response = Integer.parseInt(cl.getOptionValue(cmd_res_number)) - 1;
		}

		if (cl.hasOption(cmd_method)) {
			linkfunction = Integer.parseInt(cl.getOptionValue(cmd_method));
		}
		if (cl.hasOption(cmd_cv)) {
			cv = Integer.parseInt(cl.getOptionValue(cmd_cv));
		}
		if (cl.hasOption(cmd_seed)) {
			seed = Integer.parseInt(cl.getOptionValue(cmd_seed));
		}
		if (cl.hasOption(cmd_simu)) {
			simu = Integer.parseInt(cl.getOptionValue(cmd_simu));
		}
		if (cl.hasOption(cmd_perm)) {
			permutation = Integer.parseInt(cl.getOptionValue(cmd_perm));
		}
		if (cl.hasOption(cmd_perm_scheme)) {
			permu_scheme = true;
		}
		if (cl.hasOption(cmd_perm_fam)) {// when --pf is pronounced,
			// permu_scheme is turned on.
			permu_fam = true;
			permu_scheme = true;
		}
		if (cl.hasOption(cmd_unrelated_only)) {
			unrelated_only = true;
		}
		if (cl.hasOption(cmd_order)) {
			order = Integer.parseInt(cl.getOptionValue(cmd_order));
		}
		if (cl.hasOption(cmd_thin)) {
			thin = Double.parseDouble(cl.getOptionValue(cmd_thin));
		}
		if (cl.hasOption(cmd_slice)) {
			String[] s = cl.getOptionValue(cmd_slice).split("/");
			slice = Integer.parseInt(s[0]);
			sliceN = Integer.parseInt(s[1]);
		}
		if (cl.hasOption(cmd_missing_phenotype)) {
			missing_phenotype = cl.getOptionValue(cmd_missing_phenotype);
		}
		if (cl.hasOption(cmd_missing_genotype)) {
			missing_genotype = cl.getOptionValue(cmd_missing_genotype);
		}
		if (cl.hasOption(cmd_missing_allele)) {
			missing_allele = cl.getOptionValue(cmd_missing_allele);
		}
		if (cl.hasOption(cmd_status_shift)) {
			status_shift = -1;
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
						System.err.println("unknown covariat " + predictor_name[i]);
						System.exit(0);
					} else {
						idx.add(new Integer(j));
					}
				}
			}
			predictor = new int[idx.size()];
			for (int i = 0; i < predictor.length; i++) {
				predictor[i] = idx.get(i).intValue();
			}
		}
	}

	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}
}
