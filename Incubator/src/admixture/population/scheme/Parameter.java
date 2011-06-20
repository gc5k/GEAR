package admixture.population.scheme;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import admixture.AdmixtureConstant;

public class Parameter {

	private final String sep=",";

	private final String cmd_path = "dir";
	protected String dir = System.getProperty("user.dir") + System.getProperty("file.separator");//System.getProperty("user.dir") + System.getProperty("file.separator");
	
	private final String cmd_sd = "sd";
	protected long seed = 2011;
	
	private final String cmd_control_chr = "ctrlchr";
	protected int control_chr = 0;

	private final String cmd_null_hypothesis = "null";
	protected boolean isNullHypothesis = false;

	private final String cmd_fam_num = "fam";
	protected int[] family = new int[] { 100 };

	private final String cmd_aff_parent = "afp";
	protected int affectedParent = 0;

	private final String cmd_kid_num = "kid";
	protected int[] kid = new int[] { 2 };

	private final String cmd_aff_kid = "afk";
	protected int[] affectedKid = new int[] { 1 };

	private final String cmd_sampling_scheme = "ss";
	protected  int samplingScheme = AdmixtureConstant.FamilyExactAffected;
	
	private final String cmd_case_num = "cs";
	protected int[] cases = new int[] { 100 };

	private final String cmd_control_num = "ctrl";
	protected int[] controls = new int[] { 100 };

	private final String cmd_pop_prop = "pc";
	protected double[] popProportion = new double[] { 1, 1 };

	private final String cmd_generation = "gn";
	protected int generation = 9;

	private final String cmd_pop_prevalence = "pv";
	protected double[] popPrevalence = new double[] { 0.1, 0.2 };

	private final String cmd_geno_fun = "gf";
	protected String[] genotypeFunction = new String[] { "1111", "2121", "2222" };

	private final String cmd_geno_effect = "ge";
	protected double[] genotypeEffect = new double[] { -1.733, -1733, -1.733 };
	
	private final String cmd_linked_chr = "lc";
	protected int[] diseaseChr = new int[] { 1, 1 };

	private final String cmd_linked_locus = "ll";
	protected int[] diseaseLocus = new int[] { 0, 3 };


	private final String cmd_mu = "u";
	protected double mu = 0;
	private final String cmd_cov_effect = "ce";
	protected double covariateEffect = 0;
	private final String cmd_cov_sd = "csd";
	protected double covariateSD = 0;
	private final String cmd_esd = "rsd";
	protected double err = 1;

	private String cmd_file = "f";
	protected String[][] AIM_file = null;
	
	private String cmd_aim = "aim";
	protected int[] aim = new int[] {0, 0};
	
	private final String cmd_rep = "simu";
	protected int simulation = 10;
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
		ops.addOption(OptionBuilder.withLongOpt("family-size").withDescription("number of nuclear families").hasArg()
				.withArgName("family size").create(cmd_fam_num));
		ops.addOption(OptionBuilder.withLongOpt("kid-number").withDescription("number of kids in each family").hasArg()
				.withArgName("litter size").create(cmd_kid_num));
		ops.addOption(OptionBuilder.withLongOpt("affected-parent-number").withDescription("number of affected parents in each family").hasArg()
				.withArgName("affected parents number").create(cmd_aff_parent));
		ops.addOption(OptionBuilder.withLongOpt("affected-kid-number").withDescription("number of affected kids in each family").hasArg()
				.withArgName("affected kid number").create(cmd_aff_kid));
		ops.addOption(OptionBuilder.withLongOpt("sampling_scheme").withDescription("scheme for sampling sibs").hasArg()
				.withArgName("Sampling scheme, integer").create(cmd_sampling_scheme));
		ops.addOption(OptionBuilder.withLongOpt("case-number").withDescription("number of cases").hasArg()
				.withArgName("case number").create(cmd_case_num));
		ops.addOption(OptionBuilder.withLongOpt("control-number").withDescription("number of controls").hasArg()
				.withArgName("control number").create(cmd_control_num));
		ops.addOption(OptionBuilder.withLongOpt("population-composition").withDescription("population percent").hasArg()
				.withArgName("population percentage").create(cmd_pop_prop));
		ops.addOption(OptionBuilder.withLongOpt("generation").withDescription("generation before produce the mapping population").hasArg()
				.withArgName("successive generations for mating").create(cmd_generation));
		ops.addOption(OptionBuilder.withLongOpt("prevalence").withDescription("disease rate for each population").hasArg()
				.withArgName("disease rate").create(cmd_pop_prevalence));
		ops.addOption(OptionBuilder.withLongOpt("genotype-function").withDescription("functional genotypes").hasArg()
				.withArgName("affected genotypes").create(cmd_geno_fun));
		ops.addOption(OptionBuilder.withLongOpt("genotype-effect").withDescription("effect for functional genotypes").hasArg()
				.withArgName("functional genotype effect").create(cmd_geno_effect));
		ops.addOption(OptionBuilder.withLongOpt("disease-chr").withDescription("the chromosome of the effected genotypes").hasArg()
				.withArgName("disease chromosome").create(cmd_linked_chr));
		ops.addOption(OptionBuilder.withLongOpt("disease-locus").withDescription("the locus of the effected genotypes").hasArg()
				.withArgName("disease locus").create(cmd_linked_locus));		
		
		ops.addOption(OptionBuilder.withLongOpt("mu").withDescription("grand mean for generalized linear model").hasArg()
				.withArgName("grand mean").create(cmd_mu));
		ops.addOption(OptionBuilder.withLongOpt("cov-effect").withDescription("covariate effect for generalized linear model").hasArg()
				.withArgName("covariate effect").create(cmd_cov_effect));
		ops.addOption(OptionBuilder.withLongOpt("cov-sd").withDescription("covariate standard deviation for generalized linear model").hasArg()
				.withArgName("covariate sd").create(cmd_cov_sd));
		ops.addOption(OptionBuilder.withLongOpt("residual-sd").withDescription("covariate standard deviation for generalized linear model").hasArg()
				.withArgName("residual sd").create(cmd_esd));
		ops.addOption(OptionBuilder.withLongOpt("aim-file").withDescription("aim file").hasArgs()
				.withArgName("aim frequency file").create(cmd_file));
		ops.addOption(OptionBuilder.withLongOpt("aim-number").withDescription("aim number for each chromosome in simulation").hasArg()
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
		if(cl.hasOption(cmd_aff_parent)) {
			affectedParent = Integer.parseInt(cl.getOptionValue(cmd_aff_parent));
		}
		if(cl.hasOption(cmd_fam_num)) {
			String[] f = cl.getOptionValue(cmd_fam_num).split(sep);
			family = new int[f.length];
			for(int i =0; i < f.length; i++) family[i] = Integer.parseInt(f[i]);
		}
		if(cl.hasOption(cmd_kid_num)) {
			String[] k = cl.getOptionValue(cmd_kid_num).split(sep);
			kid = new int[k.length];
			for(int i = 0; i < k.length; i++) kid[i] = Integer.parseInt(k[i]);
		}
		if(cl.hasOption(cmd_aff_kid)) {
			String[] ak = cl.getOptionValue(cmd_aff_kid).split(sep);
			affectedKid = new int[ak.length];
			for(int i = 0; i < ak.length; i++) affectedKid[i] = Integer.parseInt(ak[i]);
		}
		if(cl.hasOption(cmd_sampling_scheme)) {
			samplingScheme = Integer.parseInt(cl.getOptionValue(cmd_sampling_scheme));
		}
		if(cl.hasOption(cmd_case_num)) {
			String[] cn = cl.getOptionValue(cmd_case_num).split(sep);
			cases = new int[cn.length];
			for(int i = 0; i < cn.length; i++) cases[i] = Integer.parseInt(cn[i]);
		}
		if(cl.hasOption(cmd_control_num)) {
			String[] cln = cl.getOptionValue(cmd_control_num).split(sep);
			controls = new int[cln.length];
			for(int i = 0; i < cln.length; i++) controls[i] = Integer.parseInt(cln[i]);
		}
		if(cl.hasOption(cmd_pop_prop)) {
			String[] pp = cl.getOptionValue(cmd_pop_prop).split(sep);
			popProportion = new double[pp.length];
			for(int i = 0; i < popProportion.length; i++) popProportion[i] = Double.parseDouble(pp[i]);
		}
		if(cl.hasOption(cmd_generation)) {
			generation = Integer.parseInt(cl.getOptionValue(cmd_generation));
		}
		if(cl.hasOption(cmd_pop_prevalence)) {
			String[] pd = cl.getOptionValue(cmd_pop_prevalence).split(sep);
			popPrevalence = new double[pd.length];
			for(int i = 0; i < pd.length; i++) popPrevalence[i] = Double.parseDouble(pd[i]);
		}
		if(cl.hasOption(cmd_geno_fun)) {
			genotypeFunction = cl.getOptionValue(cmd_geno_fun).split(sep);
		}
		if(cl.hasOption(cmd_geno_effect)) {
			String[] ge = cl.getOptionValue(cmd_geno_effect).split(sep);
			genotypeEffect = new double[ge.length];
			for(int i = 0; i <ge.length; i++) genotypeEffect[i] = Double.parseDouble(ge[i]);
		}
		if(cl.hasOption(cmd_linked_chr)) {
			String[] lc = cl.getOptionValue(cmd_linked_chr).split(sep);
			diseaseChr = new int[lc.length];
			for(int i = 0; i < lc.length; i++) diseaseChr[i] = Integer.parseInt(lc[i]);
		}
		if(cl.hasOption(cmd_linked_locus)) {
			String[] ll = cl.getOptionValue(cmd_linked_locus).split(sep);
			diseaseLocus = new int[ll.length];
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
			simulation = Integer.parseInt(cl.getOptionValue(cmd_rep));
		}
		if(cl.hasOption(cmd_file)) {
			String[] file = cl.getOptionValues(cmd_file);
			AIM_file = new String[file.length][];
			for(int i = 0, len = file.length; i < len; i++) AIM_file[i] = file[i].split(",");
		}
		if(cl.hasOption(cmd_aim)) {
			String[] aim_num = cl.getOptionValue(cmd_aim).split(sep);
			aim = new int[aim_num.length];
			for(int i = 0; i < aim_num.length; i++) aim[i] = Integer.parseInt(aim_num[i]);
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
		sb.append("Affected Parent: ");
		sb.append(affectedParent);
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
		
		sb.append("Sampling Scheme: ");
		if(samplingScheme == AdmixtureConstant.FamilyExactAffected) {
			sb.append("exact affected sibs as specified");
		} else if(samplingScheme == AdmixtureConstant.FamilyMoreThanAffected) {
			sb.append("not less than the number of affected kids as specified");
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
			for(int j = 0, len1 = AIM_file[i].length; j < len1; j++) {
				sb.append(AIM_file[i][j] + " ");
			}
			sb.append(System.getProperty("line.separator"));
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("AIM number: ");
		for(int i = 0, len = aim.length; i < len; i++) {
			sb.append(aim[i] + " ");
		}
		sb.append(System.getProperty("line.separator"));
		sb.append("Simulation Replication: " + simulation + System.getProperty("line.separator"));
		return sb.toString();
	}

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}
}
