package admixture.parameter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Parameter {

	private final String cmd_missing_allele = "missing_allele";
	public static String missing_allele = "0";

	private final String cmd_missing_phenotype = "missing_phenotype";
	public static String missing_phenotype = "99";

	private final String cmd_missing_genotype = "missing_genotype";
	public static String missing_genotype = "0";

	private final String cmd_status_shift = "1";
	public static int status_shift = 0;

	private final String cmd_mode = "md";// u, f (sii), pi(ajhg2008);
	public static String mode = "u";

	//file set start
	private final String cmd_file = "file";
	public static boolean file = false;
	private final String cmd_ped = "ped";
	public static String ped = null;
	private final String cmd_map = "map";
	public static String map = null;
	//file set end
	
	//bfile set start
	private final String cmd_bfile = "bfile";
	public static boolean bfile = false;
	private final String cmd_bed = "bed";
	public static String bed = null;
	private final String cmd_bim = "bim";
	public static String bim = null;
	private final String cmd_fam = "fam";
	public static String fam = null;
	//bfile set end

	//tfile set start
	private final String cmd_tfile = "tfile";
	public static boolean tfile = false;
	private final String cmd_tped = "tped";
	public static String tped = null;
	private final String cmd_tfam = "tfam";
	public static String tfam = null;
	//tfile set end
	
	// phenotype set start
	private final String cmd_pheno = "pheno";
	public static String pheno = null;

	private final String cmd_res = "y";
	public int response = -1;

	private final String cmd_covar_number = "covar_number";
	public static int[] predictor = null;

	private final String cmd_covar_name = "covar_name";
	public static String[] predictor_name = null;
	// phenotype set end

	// snp selection start
	private final String cmd_includesnp = "includesnp";
	public String[] includesnp = null;

	private final String cmd_includesnp_f = "includesnp_file";
	public String includesnp_file = null;

	private final String cmd_excludesnp = "excludesnp";
	public String[] excludesnp = null;

	private final String cmd_excludesnp_f = "excludesnp_file";
	public String excludesnp_file = null;

	private final String cmd_snpfile = "snpfile";
	public static String snpfile = null;

	private final String cmd_interactionsnp = "interactionsnp";
	public static String[] interactionsnp = null;

	private final String cmd_interactionsnp_f = "interactionsnp_file";
	public static String[] interactionsnp_file = null;
	// snp selection end

	// private final String cmd_topN = "topN";
	// public int topN = 5;

	private final String cmd_header = "header";
	public static boolean header = false;

	private final String cmd_method = "model";
	public int linkfunction = 0;
	private String[] lf = new String[] { "0 for linear regression", "1 for logistic regression" };

	private final String cmd_cv = "cv";
	public static int cv = 5;

	private final String cmd_order = "order";
	public int order = 1;

	// private final String cmd_min = "min";
	// public int min = 1;

	// private final String cmd_max = "max";
	// public int max = 1;

	private final String cmd_sd = "seed";
	public static int seed = 2011;

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
		ops.addOption(OptionBuilder.withDescription("specify 1 or more covariates by number.").hasArgs().create(cmd_covar_number));
		ops.addOption(OptionBuilder.withDescription("specify 1 or more covariates by name.").hasArgs().create(cmd_covar_name));

		ops.addOption(OptionBuilder.withDescription("include snps in detecting interaction").hasArgs().create(cmd_includesnp));
		ops.addOption(OptionBuilder.withDescription("specify the file containing included snps when detecting interaction").hasArg().create(
				cmd_includesnp_f));
		ops.addOption(OptionBuilder.withDescription("exclude snps in detecting interaction").hasArgs().create(cmd_excludesnp));
		ops.addOption(OptionBuilder.withDescription("specify the file containing excluded snps when detecting interaction").hasArg().create(
				cmd_excludesnp_f));

		ops.addOption(OptionBuilder.withDescription("specify interacting snps").hasArgs().create(cmd_interactionsnp));
		ops.addOption(OptionBuilder.withDescription("specify the files containing interacting snps").hasArgs().create(cmd_interactionsnp_f));

		ops.addOption(OptionBuilder.withDescription("interaction snp file").hasArg().create(cmd_snpfile));
		// ops.addOption(OptionBuilder.withDescription("only keep the results for the top N combinations").hasArg().create(cmd_topN));
		ops.addOption(OptionBuilder.withDescription("index for the response excluding the first two columns.").hasArg().create(cmd_res));

		ops.addOption(OptionBuilder.withDescription(
				"method for adjustment of the phenotype, 0 (default) for linear regression, 1 for logistic regression.").hasArg().create(cmd_method));
		ops.addOption(OptionBuilder.withDescription("fold of cross-validation, and default is 5.").hasArg().create(cmd_cv));
		ops.addOption(OptionBuilder.withDescription("specify the order of interaction.").hasArg().create(cmd_order));
		// ops.addOption(OptionBuilder.withDescription("minimal order of the interaction being searched for.").hasArg().create(cmd_min));
		// ops.addOption(OptionBuilder.withDescription("maximal order of the interaction being searched for.").hasArg().create(cmd_max));
		ops.addOption(OptionBuilder.withDescription("seed for the algorithms").hasArg().create(cmd_sd));

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
		if (cl.hasOption(cmd_mode)) {
			mode = cl.getOptionValue(cmd_mode);
		}
//file
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
		if(ped != null && map != null) {
			file = true;
		}

//tfile
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
			tfile = true;
		}

//bfile
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
			bfile = true;
		}

		if (cl.hasOption(cmd_pheno)) {
			pheno = cl.getOptionValue(cmd_pheno);
		}
		if (cl.hasOption(cmd_covar_number)) {
			String[] p = cl.getOptionValues(cmd_covar_number);
			ArrayList<Integer> idx = NewIt.newArrayList();
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
			for (int i = 0; i < predictor.length; i++)
				predictor[i] = idx.get(i).intValue();
		}
		if (cl.hasOption(cmd_covar_name)) {
			predictor_name = cl.getOptionValues(cmd_covar_name);
			for (int i = 0; i < predictor_name.length; i++) {
				if (predictor_name[i].contains("-")) {
					String[] pp = predictor_name[i].split("-");
					if (pp.length != 2) {
						System.err.println("unknown parameter " + predictor_name[i]);
						System.exit(0);
					}
				}
			}
		}
		if (cl.hasOption(cmd_includesnp)) {
			includesnp = cl.getOptionValues(cmd_includesnp);
		}
		if (cl.hasOption(cmd_includesnp_f)) {
			includesnp_file = cl.getOptionValue(cmd_includesnp_f);
			File f = new File(includesnp_file);
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(f));
			} catch (IOException E) {
				System.err.println("can't open map file.");
			}
			ArrayList<String> snplist = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					snplist.add(line);
				}
				reader.close();
			} catch (IOException E) {
				System.err.println("bad includesnp file.");
				System.exit(0);
			}
			if(snplist.size()>0) {
				includesnp = (String[]) snplist.toArray(new String[0]);
			}
		}
		if (cl.hasOption(cmd_excludesnp)) {
			excludesnp = cl.getOptionValues(cmd_excludesnp);
		}
		if (cl.hasOption(cmd_excludesnp_f)) {
			excludesnp_file = cl.getOptionValue(cmd_excludesnp_f);
			File f = new File(excludesnp_file);
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(f));
			} catch (IOException E) {
				System.err.println("can't open map file.");
			}
			ArrayList<String> snplist = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					snplist.add(line);
				}
				reader.close();
			} catch (IOException E) {
				System.err.println("bad includesnp file.");
				System.exit(0);
			}
			if(snplist.size()>0) {
				excludesnp = (String[]) snplist.toArray(new String[0]);
			}
		}

		if (cl.hasOption(cmd_header)) {
			header = true;
		}
		if (cl.hasOption(cmd_interactionsnp)) {
			interactionsnp = cl.getOptionValues(cmd_interactionsnp);
		}
		if (cl.hasOption(cmd_interactionsnp_f)) {
			interactionsnp_file = cl.getOptionValues(cmd_interactionsnp_f);
		}
		if (cl.hasOption(cmd_snpfile)) {
			snpfile = cl.getOptionValue(cmd_snpfile);
		}
		// if (cl.hasOption(cmd_topN)) {
		// topN = Integer.parseInt(cl.getOptionValue(cmd_topN));
		// }
		if (cl.hasOption(cmd_res)) {
			response = Integer.parseInt(cl.getOptionValue(cmd_res)) - 1;
		}

		if (cl.hasOption(cmd_method)) {
			linkfunction = Integer.parseInt(cl.getOptionValue(cmd_method));
		}
		if (cl.hasOption(cmd_cv)) {
			cv = Integer.parseInt(cl.getOptionValue(cmd_cv));
		}
		if (cl.hasOption(cmd_sd)) {
			seed = Integer.parseInt(cl.getOptionValue(cmd_sd));
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
		// if (cl.hasOption(cmd_min)) {
		// min = Integer.parseInt(cl.getOptionValue(cmd_min));
		// }
		// if (cl.hasOption(cmd_max)) {
		// max = Integer.parseInt(cl.getOptionValue(cmd_max));
		// }
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
					int j1 = Arrays.binarySearch(F, pp[0]);
					int j2 = Arrays.binarySearch(F, pp[1]);
					if (j1 > j2) {
						j1 = j1 ^ j2;
						j2 = j1 ^ j2;
						j1 = j1 ^ j2;
					}
					for (int j = j1; j <= j2; j++) {
						idx.add(new Integer(j));
					}
				} else {
					int j = Arrays.binarySearch(F, predictor_name[i]);
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
