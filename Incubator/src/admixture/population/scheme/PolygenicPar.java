package admixture.population.scheme;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class PolygenicPar {

	private final String sep=",";

	private final String cmd_marker = "marker";
	public static int marker = 1000;
	
	private final String cmd_LD = "ld";
	public static double ld = 0;

	private final String cmd_seed = "seed";
	public static long seed = 2011;

	private final String cmd_no_correct = "no_correct";
	private final String cmd_no_correct_long = "no-correct";
	public static boolean correct = true;

	private final String cmd_sample = "sample_size";
	private final String cmd_sample_long = "sample-size";
	public static int sample = 1000;

	private final String cmd_case = "case";
	public static int cs = 500;

	private final String cmd_prevalence = "K";
	public static double K = 0.05;

	private final String cmd_h2 = "h2";
	public static double h2 = 0.5;

	private final String cmd_out = "out";
	public static String out = "Poly";

	private final String cmd_help = "help";

	private Options ops = new Options();
	private CommandLineParser parser = new PosixParser();

	public PolygenicPar() {
		commandInitial();
	}

	public Options getOptions() {
		return ops;
	}

	public void commandInitial() {

		ops.addOption(OptionBuilder.withDescription("number of markers, defualt= " + marker).hasArg().create(cmd_marker));

		ops.addOption(OptionBuilder.withDescription("LD (correlation), defualt= " + ld).hasArg().create(cmd_LD));

		ops.addOption(OptionBuilder.withDescription("seed for simulation, default = " + seed).hasArg().withArgName("seed").create(cmd_seed));

		ops.addOption(OptionBuilder.withLongOpt(cmd_no_correct_long).withDescription("correct mean for liability, default = " + correct).hasArg().withArgName("correct").create(cmd_no_correct));

		ops.addOption(OptionBuilder.withLongOpt(cmd_sample_long).withDescription("sample size, default = " + sample).hasArg().create(cmd_sample));

		ops.addOption(OptionBuilder.withDescription("number of cases, default = " + cs).hasArg().create(cmd_case));

		ops.addOption(OptionBuilder.withDescription("prevalence, default = " + K).hasArg().create(cmd_prevalence));

		ops.addOption(OptionBuilder.withDescription("heritability, default = " + h2).hasArg().create(cmd_h2));

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
			formatter.printHelp("UGMDR", ops);
			System.exit(1);
		}
		if(cl.hasOption(cmd_marker)) {
			marker = Integer.parseInt(cl.getOptionValue(cmd_marker));
		}
		
		if(cl.hasOption(cmd_LD)) {
			ld = Double.parseDouble(cl.getOptionValue(cmd_LD));
		}

		if(cl.hasOption(cmd_seed)) {
			seed = Long.parseLong(cl.getOptionValue(cmd_seed));
		}
		
		if(cl.hasOption(cmd_no_correct)) {
			correct = false;
		}

		if(cl.hasOption(cmd_sample)) {
			sample = Integer.parseInt(cl.getOptionValue(cmd_sample)); 
		}
		
		if(cl.hasOption(cmd_case)) {
			cs = Integer.parseInt(cl.getOptionValue(cmd_case));
		}
		
		if(cl.hasOption(cmd_prevalence)) {
			K = Double.parseDouble(cl.getOptionValue(cmd_prevalence));
		}
		
		if(cl.hasOption(cmd_h2)) {
			h2 = Double.parseDouble(cl.getOptionValue(cmd_h2));
		}

		if(cl.hasOption(cmd_out)) {
			out = cl.getOptionValue(cmd_out);
		}

	}

	public static void main(String[] args) {
		PolygenicPar p = new PolygenicPar();
		p.commandListenor(args);
		System.out.println(p);
	}
}
