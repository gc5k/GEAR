package admixture.parameter;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class Parameter {
	private final String sep = ",";

	private final String cmd_mode = "md";
	public String mode = "f";

	private final String cmd_ped = "ped";
	public String pedigree = null;

	private final String cmd_phe = "phe";
	public String phenotype = null;

	private final String cmd_res = "rps";
	public int response = -1;

	private final String cmd_pred = "prd";
	public int[] predictor = null;

	private final String cmd_method = "mtd";
	public int linkfunction = 0;
	private String[] lf = new String[] {"0 for linear regression", "1 for logistic regression"};

	private final String cmd_cv = "cv";
	public int cv = 5;
	
	private final String cmd_min = "min";
	public int min = 1;
	
	private final String cmd_max = "max";
	public int max = 1;
	
	private final String cmd_sd = "sd";
	public int seed = 2011;

	private final String cmd_perm = "pm";
	public int permutation = 100;

	private final String cmd_perm_scheme = "ps";
	public boolean permu_scheme = false;

	private final String cmd_unrelated_only = "ur";
	public boolean unrelated_only = false;
	
	private final String cmd_simu = "simu";
	public int simu = 1;
	
	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();

	public Parameter() {
		commandInitial();
	}

	public void commandInitial() {
		ops.addOption(OptionBuilder.withLongOpt("mode").withDescription("u for unified framework and f for family study").hasArg().withArgName(
				"mode of the method").create(cmd_mode));
		ops.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file").hasArg().withArgName("the name of the pedigree file")
				.create(cmd_ped));
		ops.addOption(OptionBuilder.withLongOpt("phenotype").withDescription("phenotype file").hasArg().withArgName("the name of the phenotype file")
				.create(cmd_phe));
		ops.addOption(OptionBuilder.withLongOpt("response").withDescription("index for the response").hasArg().withArgName("response variable index")
				.create(cmd_res));
		ops.addOption(OptionBuilder.withLongOpt("predictor").withDescription("index(es) for the predictors").hasArg().withArgName(
				"predictor variable index(es)").create(cmd_pred));
		ops.addOption(OptionBuilder.withLongOpt("method").withDescription("method for adjustment").hasArg().withArgName(
				"generalized linear regression").create(cmd_method));
		ops.addOption(OptionBuilder.withLongOpt("cross-validation").withDescription("fold for cross-validation").hasArg().withArgName("fold for cv")
				.create(cmd_cv));
		ops.addOption(OptionBuilder.withLongOpt("minimal").withDescription("minimal order of the interaction").hasArg().withArgName("min")
				.create(cmd_min));
		ops.addOption(OptionBuilder.withLongOpt("maximal").withDescription("maximal order of the interaction").hasArg().withArgName("max")
				.create(cmd_max));
		ops.addOption(OptionBuilder.withLongOpt("seed").withDescription("seed for the algorithms").hasArg().withArgName("seed").create(cmd_sd));
		ops.addOption(OptionBuilder.withLongOpt("permutation").withDescription("replication for permutation").hasArg().withArgName("permutation")
				.create(cmd_perm));
		ops.addOption(OptionBuilder.withLongOpt("perm-scheme").withDescription("replication for permutation").create(cmd_perm_scheme));
		ops.addOption(OptionBuilder.withLongOpt("unrelated").withDescription("use unrelated indivuduals only").create(cmd_unrelated_only));
		ops.addOption(OptionBuilder.withLongOpt("simulation").withDescription("this parameter is for simulation").hasArg().withArgName(
		"replications for simulation").create(cmd_simu));
	}

	public void commandListenor(String[] args) {
		CommandLine cl = null;
		try {
			cl = parser.parse(ops, args);
		} catch (ParseException E) {
			E.printStackTrace(System.err);
		}
		if (cl.hasOption(cmd_mode)) {
			mode = cl.getOptionValue(cmd_mode);
		}
		if (cl.hasOption(cmd_ped)) {
			pedigree = cl.getOptionValue(cmd_ped);
		}
		if (cl.hasOption(cmd_phe)) {
			phenotype = cl.getOptionValue(cmd_phe);
		}
		if (cl.hasOption(cmd_res)) {
			response = Integer.parseInt(cl.getOptionValue(cmd_res)) - 1;
		}
		if (cl.hasOption(cmd_pred)) {
			String[] p = cl.getOptionValue(cmd_pred).split(sep);
			predictor = new int[p.length];
			for (int i = 0, len = p.length; i < len; i++) predictor[i] = Integer.parseInt(p[i]) - 1;
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
		if (cl.hasOption(cmd_unrelated_only)) {
			unrelated_only = true;
		}
		if (cl.hasOption(cmd_min)) {
			min = Integer.parseInt(cl.getOptionValue(cmd_min));
		}
		if (cl.hasOption(cmd_max)) {
			max = Integer.parseInt(cl.getOptionValue(cmd_max));
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("seed: ");
		sb.append(seed);
		sb.append(System.getProperty("line.separator"));

		sb.append("mode: ");
		sb.append(mode);
		sb.append(System.getProperty("line.separator"));

		sb.append("ped: ");
		sb.append(pedigree);
		sb.append(System.getProperty("line.separator"));

		sb.append("phe: ");
		sb.append(phenotype);
		sb.append(System.getProperty("line.separator"));

		sb.append("response index: ");
		sb.append(response);
		sb.append(System.getProperty("line.separator"));

		sb.append("predictor index(es): ");
		if (predictor != null) {
			for(int i = 0, len = predictor.length; i < len; i++) {
				sb.append(predictor[i] + " ");
			}
		} else {
			sb.append("none");
		}
		sb.append(System.getProperty("line.separator"));

		sb.append("link function: ");
		sb.append(lf[linkfunction]);
		sb.append(System.getProperty("line.separator"));

		sb.append("cross-validation: ");
		sb.append(cv);
		sb.append(System.getProperty("line.separator"));

		sb.append("min: ");
		sb.append(min);
		sb.append(System.getProperty("line.separator"));

		sb.append("max: ");
		sb.append(max);
		sb.append(System.getProperty("line.separator"));

		sb.append("permutation: ");
		sb.append(permutation);
		sb.append(System.getProperty("line.separator"));

		sb.append("permutation scheme: ");
		sb.append(permu_scheme);
		sb.append(System.getProperty("line.separator"));
		
		sb.append("use unrelated individuals only: ");
		sb.append(unrelated_only);
		sb.append(System.getProperty("line.separator"));

		sb.append("simulation replication: ");
		sb.append(simu);
		sb.append(System.getProperty("line.separator"));

		return sb.toString();
	}
	
	public static void main(String[] args) throws IOException {
		Parameter p = new Parameter();
		p.commandListenor(args);
		System.out.println(p);
	}
}
