package admixture.population.scheme;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class Parameter {
	private final String cmd_sd = "sd";
	private long seed = 2011;
	
	private final String cmd_control_chr = "cc";
	private int control_chr = 0;
	
	private final String cmd_null_hypothesis = "nh";
	private boolean isNullHypothesis = false;

	private final String cmd_fam_num = "fn";
	private int[] family = new int[] { 100 };
	
	private final String cmd_kid_num = "kn";
	private int[] kid = new int[] { 2 };

	private final String cmd_aff_kid = "akn";
	private int[] affectedKid = new int[] { 1 };
	
	private final String cmd_case_num = "cn";
	private int[] cases = new int[] { 100 };
	
	private final String cmd_control_num = "cln";
	private int[] controls = new int[] { 100 };
	
	private final String cmd_pop_prop = "pp";
	private double[] popProportion = new double[] { 1 };
	
	private final String cmd_pop_disease_rate = "dr";
	private double[] popDiseaseRate = new double[] { 0 };

	private final String cmd_geno_fun = "gf";
	private String[] genotypeFunction = new String[] { "0000", "1010", "1111" };
	
	private final String cmd_geno_effect = "ge";
	private double[] genotypeEffect = new double[] { 1, 1, 1 };
	
	private final String cmd_linked_chr = "lc";
	private int[] diseaseChr = new int[] { 1, 1 };

	private final String cmd_linked_locus = "ll";
	private int[] diseaseLocus = new int[] { 1, 3 };

	private final String cmd_mu = "u";
	private double mu = 0;
	private final String cmd_cov_effect = "ce";
	private double covariateEffect = 0;
	private final String cmd_cov_sd = "sd";
	private double covariateSD = 0;
	private final String cmd_esd = "rsd";
	private double err = 1;

	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();
	public Parameter() {
		commandInitial();
	}

	public void commandInitial() {
		ops.addOption(OptionBuilder.withLongOpt("seed").withDescription("seed for simulation").hasArg()
				.withArgName("seed").create(cmd_sd));
		ops.addOption(OptionBuilder.withLongOpt("control-chr").withDescription("control chromosome").hasArg()
				.withArgName("control chromosome").create(cmd_control_chr));
		ops.addOption(OptionBuilder.withLongOpt("isnullhypothesis").withDescription("generate the sample under the null hypothesis").hasArg(false)
				.withArgName("is null hypothesis").create(cmd_null_hypothesis));
		ops.addOption(OptionBuilder.withLongOpt("family-size").withDescription("number of nuclear families").hasArgs()
				.withArgName("family size").create(cmd_fam_num));
		ops.addOption(OptionBuilder.withLongOpt("kid-number").withDescription("number of kids in each family").hasArgs()
				.withArgName("litter size").create(cmd_kid_num));
		ops.addOption(OptionBuilder.withLongOpt("affected-kid-number").withDescription("number of affected kids in each family").hasArgs()
				.withArgName("affected kid number").create(cmd_aff_kid));
		ops.addOption(OptionBuilder.withLongOpt("case-number").withDescription("number of cases").hasArgs()
				.withArgName("case number").create(cmd_case_num));
		ops.addOption(OptionBuilder.withLongOpt("control-number").withDescription("number of controls").hasArgs()
				.withArgName("control number").create(cmd_control_num));
		ops.addOption(OptionBuilder.withLongOpt("population-composition").withDescription("population percent").hasArgs()
				.withArgName("population percentage").create(cmd_pop_prop));
		ops.addOption(OptionBuilder.withLongOpt("disease-rate").withDescription("disease rate for each population").hasArgs()
				.withArgName("disease rate").create(cmd_pop_disease_rate));
		ops.addOption(OptionBuilder.withLongOpt("genotype-function").withDescription("functional genotypes").hasArgs()
				.withArgName("affected genotypes").create(cmd_geno_fun));
		ops.addOption(OptionBuilder.withLongOpt("genotype-effect").withDescription("effect for functional genotypes").hasArgs()
				.withArgName("functional genotype effect").create(cmd_geno_effect));
		ops.addOption(OptionBuilder.withLongOpt("disease-chr").withDescription("the chromosome of the effected genotypes").hasArgs()
				.withArgName("disease chromosome").create(cmd_linked_chr));
		ops.addOption(OptionBuilder.withLongOpt("disease-locus").withDescription("the locus of the effected genotypes").hasArgs()
				.withArgName("disease locus").create(cmd_linked_locus));		
		
		ops.addOption(OptionBuilder.withLongOpt("mu").withDescription("grand mean for generalized linear model").hasArg()
				.withArgName("grand mean").create(cmd_mu));
		ops.addOption(OptionBuilder.withLongOpt("cov-effect").withDescription("covariate effect for generalized linear model").hasArg()
				.withArgName("covariate effect").create(cmd_cov_effect));
		ops.addOption(OptionBuilder.withLongOpt("cov-sd").withDescription("covariate standard deviation for generalized linear model").hasArg()
				.withArgName("covariate sd").create(cmd_cov_sd));
		ops.addOption(OptionBuilder.withLongOpt("residual-sd").withDescription("covariate standard deviation for generalized linear model").hasArg()
				.withArgName("residual sd").create(cmd_cov_sd));
	}

	public void commandListenor(String[] args) {
		CommandLine cl = null;
		try {
			cl = parser.parse(ops, args);
		} catch (ParseException E) {
			E.printStackTrace(System.err);
		}
		if(cl.hasOption(cmd_sd)) {
			seed = Long.parseLong(cl.getOptionValue(cmd_sd));
		}
		if(cl.hasOption(cmd_control_chr)) {
			control_chr = Integer.parseInt(cl.getOptionValue(cmd_control_chr));
		}
		if(cl.hasOption(cmd_null_hypothesis)) {
			isNullHypothesis = true;
		}
		if(cl.hasOption(cmd_fam_num)) {
			String[] f = cl.getOptionValues(cmd_fam_num);
			family = new int[f.length];
			for(int i =0; i < f.length; i++) family[i] = Integer.parseInt(f[i]);
		}
		if(cl.hasOption(cmd_kid_num)) {
			String[] k = cl.getOptionValues(cmd_kid_num);
			kid = new int[k.length];
			for(int i = 0; i < k.length; i++) kid[i] = Integer.parseInt(k[i]);
		}
		if(cl.hasOption(cmd_aff_kid)) {
			String[] ak = cl.getOptionValues(cmd_aff_kid);
			affectedKid = new int[ak.length];
			for(int i = 0; i < ak.length; i++) affectedKid[i] = Integer.parseInt(ak[i]);
		}
		if(cl.hasOption(cmd_case_num)) {
			String[] cn = cl.getOptionValues(cmd_case_num);
			cases = new int[cn.length];
			for(int i = 0; i < cn.length; i++) cases[i] = Integer.parseInt(cn[i]);
		}
		if(cl.hasOption(cmd_control_num)) {
			String[] cln = cl.getOptionValues(cmd_control_num);
			controls = new int[cln.length];
			for(int i = 0; i < cln.length; i++) controls[i] = Integer.parseInt(cln[i]);
		}
		if(cl.hasOption(cmd_pop_prop)) {
			String[] pp = cl.getOptionValues(cmd_pop_prop);
			popProportion = new double[pp.length];
			for(int i = 0; i < popProportion.length; i++) popProportion[i] = Double.parseDouble(pp[i]);
		}
		if(cl.hasOption(cmd_pop_disease_rate)) {
			String[] pd = cl.getOptionValues(cmd_pop_disease_rate);
			popDiseaseRate = new double[pd.length];
			for(int i = 0; i < pd.length; i++) popDiseaseRate[i] = Double.parseDouble(pd[i]);
		}
		if(cl.hasOption(cmd_geno_fun)) {
			String[] genotypeFunction = cl.getOptionValues(cmd_geno_fun);
		}
		if(cl.hasOption(cmd_geno_effect)) {
			String[] ge = cl.getOptionValues(cmd_geno_fun);
			for(int i = 0; i <ge.length; i++) genotypeEffect[i] = Double.parseDouble(ge[i]);
		}
		if(cl.hasOption(cmd_linked_chr)) {
			String[] lc = cl.getOptionValues(cmd_linked_chr);
			for(int i = 0; i < lc.length; i++) diseaseChr[i] = Integer.parseInt(lc[i]);
		}
		if(cl.hasOption(cmd_linked_locus)) {
			String[] ll = cl.getOptionValues(cmd_linked_locus);
			for(int i = 0; i < ll.length; i++) diseaseLocus[i] = Integer.parseInt(ll[i]);
		}
		if(cl.hasOption(cmd_mu)) {
			mu = Double.parseDouble(cl.getOptionValue(cmd_mu));
		}
		if(cl.hasOption(cmd_cov_effect)) {
			covariateEffect = Double.parseDouble(cl.getOptionValue(cmd_cov_effect));
		}
		if(cl.hasOption(cmd_cov_sd)) {
			covariateSD = Double.parseDouble(cl.getOptionValue(cmd_cov_sd));
		}
		if(cl.hasOption(cmd_esd)) {
			err = Double.parseDouble(cl.getOptionValue(cmd_esd));
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("seed=" + seed + "\n");
		sb.append("control_chr=" + control_chr + "\n");
		sb.append("is null hypothesis: " + isNullHypothesis + "\n");
		return sb.toString();
	}

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}
}
