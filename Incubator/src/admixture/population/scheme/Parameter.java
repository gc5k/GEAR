package admixture.population.scheme;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class Parameter {

	private final String cmd_path = "dir";
	protected String dir = System.getProperty("user.dir") + System.getProperty("file.separator");//System.getProperty("user.dir") + System.getProperty("file.separator");
	
	private final String cmd_sd = "sd";
	protected long seed = 2011;
	
	private final String cmd_control_chr = "ctrlchr";
	protected int control_chr = 0;

	private final String cmd_null_hypothesis = "null";
	protected boolean isNullHypothesis = false;

	private final String cmd_fam_num = "famnum";
	protected int[] family = new int[] { 100 };

	private final String cmd_kid_num = "kidnum";
	protected int[] kid = new int[] { 2 };

	private final String cmd_aff_kid = "afkidnum";
	protected int[] affectedKid = new int[] { 1 };

	private final String cmd_case_num = "cases";
	protected int[] cases = new int[] { 100 };

	private final String cmd_control_num = "ctrl";
	protected int[] controls = new int[] { 100 };

	private final String cmd_pop_prop = "popprop";
	protected double[] popProportion = new double[] { 0.975, 0.25 };

	private final String cmd_generation = "gr";
	protected int generation = 9;

	private final String cmd_pop_prevalence = "prv";
	protected double[] popPrevalence = new double[] { 0.1, 0.2 };

	private final String cmd_geno_fun = "gf";
	protected String[] genotypeFunction = new String[] { "0000", "1010", "1111" };

	private final String cmd_geno_effect = "ge";
	protected double[] genotypeEffect = new double[] { 1, 1, 1 };
	
	private final String cmd_linked_chr = "lc";
	protected int[] diseaseChr = new int[] { 1, 1 };

	private final String cmd_linked_locus = "ll";
	protected int[] diseaseLocus = new int[] { 1, 3 };


	private final String cmd_mu = "u";
	protected double mu = 0;
	private final String cmd_cov_effect = "ce";
	protected double covariateEffect = 0;
	private final String cmd_cov_sd = "csd";
	protected double covariateSD = 0;
	private final String cmd_esd = "rsd";
	protected double err = 1;

	private String cmd_file = "f";
	protected String[][] AIM_file = new String[][]{ { "allele_freq_chr1_200snp.txt", "allele_freq_chr2.txt" }};
	
	private String cmd_aim = "aim";
	protected int[] aim = new int[] {0, 0};
	
	private final String cmd_rep = "rep";
	protected int replication = 10;
	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();

	public Parameter() {
		commandInitial();
	}

	public void commandInitial() {
		ops.addOption(OptionBuilder.withLongOpt("working-dir").withDescription("working directory").hasArg()
				.withArgName("working dir").create(cmd_path));		
		ops.addOption(OptionBuilder.withLongOpt("seed").withDescription("seed for simulation").hasArg()
				.withArgName("seed").create(cmd_sd));
		ops.addOption(OptionBuilder.withLongOpt("control-chr").withDescription("control chromosome").hasArg()
				.withArgName("control chromosome").create(cmd_control_chr));
		ops.addOption(OptionBuilder.withLongOpt("nullhypothesis").withDescription("generate the sample under the null hypothesis").hasArg(false)
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
		ops.addOption(OptionBuilder.withLongOpt("generation").withDescription("generation before produce the mapping population").hasArgs()
				.withArgName("successive generations for mating").create(cmd_generation));
		ops.addOption(OptionBuilder.withLongOpt("disease-rate").withDescription("disease rate for each population").hasArgs()
				.withArgName("disease rate").create(cmd_pop_prevalence));
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
		ops.addOption(OptionBuilder.withLongOpt("aim-file").withDescription("aim file").hasArgs()
				.withArgName("aim frequency file").create(cmd_file));
		ops.addOption(OptionBuilder.withLongOpt("aim-number").withDescription("aim number for each chromosome in simulation").hasArgs()
				.withArgName("number of aim markers").create(cmd_aim));
		
		ops.addOption(OptionBuilder.withLongOpt("replication").withDescription("replications for simulation").hasArg()
				.withArgName("replications").create(cmd_rep));
	}

	
	public void commandListenor(String[] args) {
		CommandLine cl = null;
		try {
			cl = parser.parse(ops, args);
		} catch (ParseException E) {
			E.printStackTrace(System.err);
		}
		if(cl.hasOption(cmd_path)) {
			dir = cl.getOptionValue(cmd_path);
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
		if(cl.hasOption(cmd_generation)) {
			generation = Integer.parseInt(cl.getOptionValue(cmd_generation));
		}
		if(cl.hasOption(cmd_pop_prevalence)) {
			String[] pd = cl.getOptionValues(cmd_pop_prevalence);
			popPrevalence = new double[pd.length];
			for(int i = 0; i < pd.length; i++) popPrevalence[i] = Double.parseDouble(pd[i]);
		}
		if(cl.hasOption(cmd_geno_fun)) {
			genotypeFunction = cl.getOptionValues(cmd_geno_fun);
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
		if(cl.hasOption(cmd_rep)) {
			replication = Integer.parseInt(cl.getOptionValue(cmd_rep));
		}
		if(cl.hasOption(cmd_file)) {
			String[] file = cl.getArgs();
			AIM_file = new String[file.length][];
			for(int i = 0, len = file.length; i < len; i++) {
				AIM_file[i] = file[i].split(",");
			}
		}
		if(cl.hasOption(cmd_aim)) {
			String[] aim_num = cl.getOptionValues(cmd_aim);
			for(int i = 0; i < aim_num.length; i++) {
				aim[i] = Integer.parseInt(aim_num[i]);
			}
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("working dir=" + dir + System.getProperty("line.separator"));
		sb.append("seed=" + seed + System.getProperty("line.separator"));
		sb.append("control_chr=" + control_chr + System.getProperty("line.separator"));
		sb.append("is null hypothesis: " + isNullHypothesis + System.getProperty("line.separator"));

		sb.append("Family: ");
		for(int i = 0, len = family.length; i<len; i++) {
			sb.append(family[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Kid: ");
		for(int i = 0, len = kid.length; i < len; i++) {
			sb.append(kid[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Affected Kid: ");
		for(int i = 0, len = affectedKid.length; i < len; i++) {
			sb.append(affectedKid[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));
		
		sb.append("Cases: ");
		for(int i = 0, len = cases.length; i < len; i++) {
			sb.append(cases[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));
		
		sb.append("Controls: ");
		for(int i = 0, len = controls.length; i < len; i++) {
			sb.append(controls[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));
		
		sb.append("Genotype Function: ");
		for(int i = 0, len = genotypeFunction.length; i < len; i++) {
			sb.append(genotypeFunction[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Genotype Affect: ");
		for(int i = 0, len = genotypeEffect.length; i < len; i++) {
			sb.append(genotypeEffect[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Disease Chromosome: ");
		for(int i = 0, len = diseaseChr.length; i < len; i++) {
			sb.append(diseaseChr[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Disease Locus: ");
		for(int i = 0, len = diseaseLocus.length; i < len; i++) {
			sb.append(diseaseLocus[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Mu: " + mu + System.getProperty("line.separator"));
		sb.append("covariate effect: " + covariateEffect + System.getProperty("line.separator"));
		sb.append("covariate sd: " + covariateSD + System.getProperty("line.separator"));
		sb.append("residual sd: " + err + System.getProperty("line.separator"));

		sb.append("Population Proportion: ");
		for(int i = 0, len = popProportion.length; i < len; i++) {
			sb.append(popProportion[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Generation: " + generation + System.getProperty("line.separator"));

		sb.append("Population Disease Rate: ");
		for(int i = 0, len = popPrevalence.length; i < len; i++) {
			sb.append(popPrevalence[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("AIM file: ");
		for(int i = 0, len = AIM_file.length; i < len; i++) {
			sb.append(AIM_file[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("Replication: " + replication + System.getProperty("line.separator"));
		return sb.toString();
	}

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
		
		System.setProperty("java.io.tempdir", "c:/cgb/");
		System.out.println(System.getProperty("user.name"));
	}
}
