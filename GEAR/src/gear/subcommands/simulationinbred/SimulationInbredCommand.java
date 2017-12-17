package gear.subcommands.simulationinbred;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SimulationInbredCommand extends Command
{

	@Override
	public String getName()
	{
		return "simuinbred";
	}

	@Override
	public String getDescription()
	{
		return "Simulation polygenic model for inbred lines";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SAMPLE_SIZE_LONG_DESC).withLongOpt(OPT_SAMPLE_SIZE_LONG).hasArg().isRequired().create(OPT_SAMPLE_SIZE));

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
		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));

		options.addOption(OptionBuilder.withDescription(OPT_RAND_LD_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_LD_RANGE_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_LD_FILE_LONG_DESC).withLongOpt(OPT_LD_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));

		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArg().create(OPT_HSQ));

		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_LONG_DESC).withLongOpt(OPT_MAKE_BED_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		
		options.addOption(OptionBuilder.withDescription(OPT_FAM_ID_PREFIX_LONG_DESC).withLongOpt(OPT_FAM_ID_PREFIX_LONG).hasArg().create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		SimulationInbredCommandArguments simuInbredArgs = new SimulationInbredCommandArguments();

		if (cmdLine.hasOption(OPT_SAMPLE_SIZE_LONG))
		{
			simuInbredArgs.setSampleSize(cmdLine.getOptionValue(OPT_SAMPLE_SIZE_LONG));
		}
		
		if (cmdLine.hasOption(OPT_REP))
		{
			simuInbredArgs.setRep(cmdLine.getOptionValue(OPT_REP));
		}

		if (cmdLine.hasOption(OPT_MARKER_LONG))
		{
			simuInbredArgs.setMarkerNum(cmdLine.getOptionValue(OPT_MARKER_LONG));
		}

		if (cmdLine.hasOption(OPT_NULL_MARKER_LONG))
		{
			simuInbredArgs.setNullMarkerNum(cmdLine.getOptionValue(OPT_NULL_MARKER_LONG));
		}

		//add eff
		simuInbredArgs.setPolyEffect();

		if (cmdLine.hasOption(OPT_EFFECT))
		{
			simuInbredArgs.setPlainEffect(parseDoubleOptionValue(cmdLine, OPT_EFFECT, "0.5"));
		}

		if (cmdLine.hasOption(OPT_POLY_EFFECT_SORT_LONG))
		{
			simuInbredArgs.setPolyEffectSort();
		}

		if (cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			simuInbredArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		//dom eff

		simuInbredArgs.setFreq(parseDoubleOptionValue(cmdLine, OPT_FREQ, "0.5"));
		
		if (cmdLine.hasOption(OPT_UNIF_FREQ_LONG))
		{
			simuInbredArgs.setUnifFreq();
		}

		if (cmdLine.hasOption(OPT_FREQ_RANGE_LONG))
		{
			simuInbredArgs.setFreqRange(cmdLine.getOptionValues(OPT_FREQ_RANGE_LONG));
		}

		if (cmdLine.hasOption(OPT_FREQ_FILE_LONG))
		{
			simuInbredArgs.setFreqFile(cmdLine.getOptionValue(OPT_FREQ_FILE_LONG));
		}
		
		simuInbredArgs.setLD(parseDoubleOptionValue(cmdLine, OPT_LD, "0.0"));

		if (cmdLine.hasOption(OPT_RAND_LD_LONG))
		{
			simuInbredArgs.setRandLD();
		}

		if (cmdLine.hasOption(OPT_LD_RANGE_LONG))
		{
			simuInbredArgs.setLDRange(cmdLine.getOptionValues(OPT_LD_RANGE_LONG));
		}

		if (cmdLine.hasOption(OPT_LD_FILE_LONG))
		{
			simuInbredArgs.setLDFile(cmdLine.getOptionValue(OPT_LD_FILE_LONG));
		}

		simuInbredArgs.setHsq(parseDoubleOptionValue(cmdLine, OPT_HSQ, "0.5"));

		if (cmdLine.hasOption(OPT_MAKE_BED_LONG))
		{
			simuInbredArgs.setMakeBed();
		}

		simuInbredArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, "2014"));
		
		if (cmdLine.hasOption(OPT_FAM_ID_PREFIX_LONG))
		{
			simuInbredArgs.setFamIDPrefix(cmdLine.getOptionValue(OPT_FAM_ID_PREFIX_LONG));			
		}
		return simuInbredArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimulationInbredCommandImpl();
	}

	private static final String OPT_SAMPLE_SIZE = "n";
	private static final String OPT_SAMPLE_SIZE_LONG = "sample-size";
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
	private static final String OPT_UNIF_FREQ_LONG_DESC = "Uniform distriubted frequency between 0.01~0.5.";

	private static final String OPT_FREQ_RANGE_LONG = "freq-range";
	private static final String OPT_FREQ_RANGE_LONG_DESC = "Uniform distriubted frequency between a~b.";

	private static final String OPT_FREQ_FILE_LONG = "freq-file";
	private static final String OPT_FREQ_FILE_LONG_DESC = "Read frequency from the file specified.";

	//ld
	private static final String OPT_LD = "ld";
	private static final String OPT_LD_DESC = "LD (Lewontin's D prime), zero by default.";

	private static final String OPT_RAND_LD_LONG = "rand-ld";
	private static final String OPT_RAND_LD_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between 0~1."; 

	private static final String OPT_LD_RANGE_LONG = "ld-range";
	private static final String OPT_LD_RANGE_LONG_DESC = "LD (Lewontin's D prime) follows uniform distribution between a~b.";

	private static final String OPT_LD_FILE_LONG = "ld-file";
	private static final String OPT_LD_FILE_LONG_DESC = "LD (Lewontin's D prime) file."; 

	//hsq
	private static final String OPT_HSQ = "hsq";
	private static final String OPT_HSQ_DESC = "Heritability for additive.";

	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_LONG_DESC = "make-bed";

	private static final String OPT_FAM_ID_PREFIX_LONG = "fam-prefix";
	private static final String OPT_FAM_ID_PREFIX_LONG_DESC = "fam-prefix.";

	private static final String OPT_REP = "rep";
	private static final String OPT_REP_DESC = "Replication for simulation. New error will be added to the same breeding values.";
}
