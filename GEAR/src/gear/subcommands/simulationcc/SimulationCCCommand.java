package gear.subcommands.simulationcc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SimulationCCCommand extends Command
{

	@Override
	public String getName()
	{
		return "simu-cc";
	}

	@Override
	public String getDescription()
	{
		return "Simulating case-control data";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SAMPLE_SIZE_LONG_DESC).withLongOpt(OPT_SAMPLE_SIZE_LONG).hasArgs(2).isRequired().create(OPT_SAMPLE_SIZE));
		options.addOption(OptionBuilder.withDescription(OPT_K_DESC).hasArg().create(OPT_K));

		options.addOption(OptionBuilder.withDescription(OPT_MARKER_LONG_DESC).withLongOpt(OPT_MARKER_LONG).hasArg().isRequired().create(OPT_MARKER));
		options.addOption(OptionBuilder.withDescription(OPT_NULL_MARKER_LONG_DESC).withLongOpt(OPT_NULL_MARKER_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_DESC).hasArg().create(OPT_EFFECT));
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_SORT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_SORT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_FREQ_DESC).hasArg().create(OPT_FREQ));
		options.addOption(OptionBuilder.withDescription(OPT_UNIF_FREQ_LONG_DESC).withLongOpt(OPT_UNIF_FREQ_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_FREQ_RANGE_LONG_DESC).withLongOpt(OPT_FREQ_RANGE_LONG).hasArgs(2).create());
		options.addOption(OptionBuilder.withDescription(OPT_FREQ_FILE_LONG_DESC).withLongOpt(OPT_FREQ_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_LD_DESC).hasArg().create(OPT_LD));

		options.addOption(OptionBuilder.withDescription(OPT_RAND_LD_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_LD_RANGE_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArg().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_LONG_DESC).withLongOpt(OPT_MAKE_BED_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		SimulationCCCommandArguments simuCCArgs = new SimulationCCCommandArguments();

		if(cmdLine.hasOption(OPT_SAMPLE_SIZE_LONG))
		{
			simuCCArgs.setSampleSize(cmdLine.getOptionValues(OPT_SAMPLE_SIZE_LONG));
		}

		if(cmdLine.hasOption(OPT_MARKER_LONG))
		{
			simuCCArgs.setMarkerNum(cmdLine.getOptionValue(OPT_MARKER_LONG));
		}

		if(cmdLine.hasOption(OPT_NULL_MARKER_LONG))
		{
			simuCCArgs.setNullMarkerNum(cmdLine.getOptionValue(OPT_NULL_MARKER_LONG));
		}

		simuCCArgs.setPolyEffect();

		if(cmdLine.hasOption(OPT_EFFECT))
		{
			simuCCArgs.setPolyEffect(parseDoubleOptionValue(cmdLine, OPT_EFFECT, "0.5"));
		}

		if(cmdLine.hasOption(OPT_POLY_EFFECT_SORT_LONG))
		{
			simuCCArgs.setPolyEffectSort();
		}

		if(cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			simuCCArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		simuCCArgs.setFreq(parseDoubleOptionValue(cmdLine, OPT_FREQ, "0.5"));
		
		if(cmdLine.hasOption(OPT_UNIF_FREQ_LONG))
		{
			simuCCArgs.setUnifFreq();
		}

		if(cmdLine.hasOption(OPT_FREQ_RANGE_LONG))
		{
			simuCCArgs.setFreqRange(cmdLine.getOptionValues(OPT_FREQ_RANGE_LONG));
		}

		if(cmdLine.hasOption(OPT_FREQ_FILE_LONG))
		{
			simuCCArgs.setFreqFile(cmdLine.getOptionValue(OPT_FREQ_FILE_LONG));
		}
		
		simuCCArgs.setLD(parseDoubleOptionValue(cmdLine, OPT_LD, "0.0"));

		if(cmdLine.hasOption(OPT_RAND_LD_LONG))
		{
			simuCCArgs.setRandLD();
		}

		if(cmdLine.hasOption(OPT_LD_RANGE_LONG))
		{
			simuCCArgs.setLDRange(cmdLine.getOptionValues(OPT_LD_RANGE_LONG));
		}

		simuCCArgs.setHsq(parseDoubleOptionValue(cmdLine, OPT_HSQ, "0.5"));

		
		if(cmdLine.hasOption(OPT_MAKE_BED_LONG))
		{
			simuCCArgs.setMakeBed();
		}
		
		if(cmdLine.hasOption(OPT_K))
		{
			simuCCArgs.setK(cmdLine.getOptionValue(OPT_K));
		}

		simuCCArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, "2014"));
		return simuCCArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimulationCCCommandImpl();
	}

	private static final String OPT_SAMPLE_SIZE = "cc";
	private static final String OPT_SAMPLE_SIZE_LONG = "case-control";
	private static final String OPT_SAMPLE_SIZE_LONG_DESC = "Sample size";

	private static final String OPT_MARKER = "m";
	private static final String OPT_MARKER_LONG = "marker";
	private static final String OPT_MARKER_LONG_DESC = "Number of markers.";

	private static final String OPT_NULL_MARKER_LONG = "null-marker";
	private static final String OPT_NULL_MARKER_LONG_DESC = "Number of null markers.";

	//effect
	private static final String OPT_EFFECT = "effect";
	private static final String OPT_EFFECT_DESC = "Equal effects, 1 by default.";

	private static final String OPT_POLY_EFFECT_LONG = "poly-effect";
	private static final String OPT_POLY_EFFECT_LONG_DESC = "Polygenic (normal distribution) model.";

	private static final String OPT_POLY_EFFECT_SORT_LONG = "poly-effect-sort";
	private static final String OPT_POLY_EFFECT_SORT_LONG_DESC = "Sorted polygenic (normal distribution) model.";

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	//freq
	private static final String OPT_FREQ = "freq";
	private static final String OPT_FREQ_DESC = "Equal frequency, 0.5 by default.";

	private static final String OPT_UNIF_FREQ_LONG = "unif-freq";
	private static final String OPT_UNIF_FREQ_LONG_DESC = "Uniform distriubted frequency between 0.01~0.5";

	private static final String OPT_FREQ_RANGE_LONG = "freq-range";
	private static final String OPT_FREQ_RANGE_LONG_DESC = "Uniform distriubted frequency between a~b";

	private static final String OPT_FREQ_FILE_LONG = "freq-file";
	private static final String OPT_FREQ_FILE_LONG_DESC = "Read frequency from the file specified.";

	//ld
	private static final String OPT_LD = "ld";
	private static final String OPT_LD_DESC = "LD (Lewontin's D prime), zero by default.";

	private static final String OPT_RAND_LD_LONG = "rand-ld";
	private static final String OPT_RAND_LD_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between 0~1"; 

	private static final String OPT_LD_RANGE_LONG = "ld-range";
	private static final String OPT_LD_RANGE_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between a~b"; 

	//hsq
	private static final String OPT_HSQ = "hsq";
	private static final String OPT_HSQ_DESC = "Heritability.";

	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_LONG_DESC = "make-bed";
	
	private static final String OPT_K = "k";
	private static final String OPT_K_DESC = "prevalence";
}
