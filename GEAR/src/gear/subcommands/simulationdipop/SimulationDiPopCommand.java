package gear.subcommands.simulationdipop;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SimulationDiPopCommand extends Command {

	@Override
	public String getName() {
		return "dipop";
	}

	@Override
	public String getDescription() {
		return "Simulating divergent simulation";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_SAMPLE_SIZE_LONG_DESC).withLongOpt(OPT_SAMPLE_SIZE_LONG)
				.hasArgs(2).isRequired().create(OPT_SAMPLE_SIZE));

		options.addOption(OptionBuilder.withDescription(OPT_MARKER_LONG_DESC).withLongOpt(OPT_MARKER_LONG).hasArgs()
				.isRequired().create(OPT_MARKER));

		options.addOption(OptionBuilder.withDescription(OPT_FREQ_RANGE_LONG_DESC).withLongOpt(OPT_FREQ_RANGE_LONG).hasArgs(2).create());

		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_DESC).hasArg().create(OPT_EFFECT));
		options.addOption(
				OptionBuilder.withDescription(OPT_POLY_EFFECT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_SORT_LONG_DESC)
				.withLongOpt(OPT_POLY_EFFECT_SORT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG)
				.hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_LD_DESC).hasArg().create(OPT_LD));

		options.addOption(OptionBuilder.withDescription(OPT_RAND_LD_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_LD_RANGE_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());

		options.addOption(
				OptionBuilder.withDescription(OPT_MAKE_BED_LONG_DESC).withLongOpt(OPT_MAKE_BED_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_FST_DESC).hasArgs().create(OPT_FST));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		SimulationDiPopCommandArguments dpArgs = new SimulationDiPopCommandArguments();

		if (cmdLine.hasOption(OPT_SAMPLE_SIZE_LONG)) {
			dpArgs.setSampleSize(cmdLine.getOptionValues(OPT_SAMPLE_SIZE_LONG));
		}

		if (cmdLine.hasOption(OPT_MARKER_LONG)) {
			dpArgs.setMarkerNum(cmdLine.getOptionValues(OPT_MARKER_LONG));
		}

		if (cmdLine.hasOption(OPT_FREQ_RANGE_LONG)) {
			dpArgs.setFreqRange(cmdLine.getOptionValues(OPT_FREQ_RANGE_LONG));
		}

		dpArgs.setPolyEffect();

		if (cmdLine.hasOption(OPT_EFFECT)) {
			dpArgs.setPlainEffect(parseDoubleOptionValue(cmdLine, OPT_EFFECT, "0.5"));
		}

		if (cmdLine.hasOption(OPT_POLY_EFFECT_SORT_LONG)) {
			dpArgs.setPolyEffectSort();
		}

		if (cmdLine.hasOption(OPT_EFFECT_FILE_LONG)) {
			dpArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		dpArgs.setLD(parseDoubleOptionValue(cmdLine, OPT_LD, "0.0"));

		if (cmdLine.hasOption(OPT_RAND_LD_LONG)) {
			dpArgs.setRandLD();
		}

		if (cmdLine.hasOption(OPT_LD_RANGE_LONG)) {
			dpArgs.setLDRange(cmdLine.getOptionValues(OPT_LD_RANGE_LONG));
		}

		dpArgs.setHsq(parseDoubleOptionValue(cmdLine, OPT_HSQ, "0.5"));

		if (cmdLine.hasOption(OPT_MAKE_BED_LONG)) {
			dpArgs.setMakeBed();
		}

		dpArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, "2014"));
		dpArgs.setFst(cmdLine.getOptionValues(OPT_FST));
		return dpArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		// TODO Auto-generated method stub
		return new SimulationDiPopCommandImpl();
	}

	private static final String OPT_SAMPLE_SIZE = "s";
	private static final String OPT_SAMPLE_SIZE_LONG = "di-size";
	private static final String OPT_SAMPLE_SIZE_LONG_DESC = "Sample size";

	private static final String OPT_MARKER = "m";
	private static final String OPT_MARKER_LONG = "marker";
	private static final String OPT_MARKER_LONG_DESC = "Number of markers.";

	// effect
	private static final String OPT_EFFECT = "effect";
	private static final String OPT_EFFECT_DESC = "Equal effects, 1 by default.";

	private static final String OPT_FREQ_RANGE_LONG = "freq-range";
	private static final String OPT_FREQ_RANGE_LONG_DESC = "Uniform distriubted frequency between a~b";

	private static final String OPT_POLY_EFFECT_LONG = "poly-effect";
	private static final String OPT_POLY_EFFECT_LONG_DESC = "Polygenic (normal distribution) model.";

	private static final String OPT_POLY_EFFECT_SORT_LONG = "poly-effect-sort";
	private static final String OPT_POLY_EFFECT_SORT_LONG_DESC = "Sorted polygenic (normal distribution) model.";

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	// ld
	private static final String OPT_LD = "ld";
	private static final String OPT_LD_DESC = "LD (Lewontin's D prime), zero by default.";

	private static final String OPT_RAND_LD_LONG = "rand-ld";
	private static final String OPT_RAND_LD_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between 0~1";

	private static final String OPT_LD_RANGE_LONG = "ld-range";
	private static final String OPT_LD_RANGE_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between a~b";

	// hsq
	private static final String OPT_HSQ = "hsq";

	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_LONG_DESC = "make-bed";

	private static final String OPT_FST = "fst";
	private static final String OPT_FST_DESC = "fst";
			
}
