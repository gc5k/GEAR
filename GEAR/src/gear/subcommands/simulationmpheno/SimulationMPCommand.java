package gear.subcommands.simulationmpheno;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SimulationMPCommand extends Command {

	@Override
	public String getName()
	{
		return "simump";
	}

	@Override
	public String getDescription()
	{
		return "simulation for correlated traits";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SAMPLE_SIZE_LONG_DESC).withLongOpt(OPT_SAMPLE_SIZE_LONG).hasArg().isRequired().create(OPT_SAMPLE_SIZE));
		options.addOption(OptionBuilder.withDescription(OPT_MARKER_LONG_DESC).withLongOpt(OPT_MARKER_LONG).hasArg().isRequired().create(OPT_MARKER));

		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_DESC).hasArg().create(OPT_EFFECT));
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_LONG).create());
//		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_SORT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_SORT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_FREQ_DESC).hasArg().create(OPT_FREQ));
		options.addOption(OptionBuilder.withDescription(OPT_UNIF_FREQ_LONG_DESC).withLongOpt(OPT_UNIF_FREQ_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_FREQ_RANGE_LONG_DESC).withLongOpt(OPT_FREQ_RANGE_LONG).hasArgs(2).create());
		options.addOption(OptionBuilder.withDescription(OPT_FREQ_FILE_LONG_DESC).withLongOpt(OPT_FREQ_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_LD_DESC).hasArg().create(OPT_LD));
		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));

		options.addOption(OptionBuilder.withDescription(OPT_RAND_LD_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_LD_RANGE_LONG_DESC).withLongOpt(OPT_RAND_LD_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArgs().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_COR_MAT_LONG_DESC).withLongOpt(OPT_COR_MAT_LONG).hasArg().create(OPT_CM));
		options.addOption(OptionBuilder.withDescription(OPT_COR_MAT_E_LONG_DESC).withLongOpt(OPT_COR_MAT_E_LONG).hasArg().create(OPT_CME));

		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_LONG_DESC).withLongOpt(OPT_MAKE_BED_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SimulationMPCommandArguments simuMPArgs = new SimulationMPCommandArguments();

		if(cmdLine.hasOption(OPT_SAMPLE_SIZE_LONG))
		{
			simuMPArgs.setSampleSize(cmdLine.getOptionValue(OPT_SAMPLE_SIZE_LONG));
		}
		
		if(cmdLine.hasOption(OPT_REP))
		{
			simuMPArgs.setRep(cmdLine.getOptionValue(OPT_REP));
		}

		if(cmdLine.hasOption(OPT_MARKER_LONG))
		{
			simuMPArgs.setMarkerNum(cmdLine.getOptionValue(OPT_MARKER_LONG));
		}

		simuMPArgs.setPolyEffect();

		if(cmdLine.hasOption(OPT_EFFECT))
		{
			simuMPArgs.setPlainEffect(parseDoubleOptionValue(cmdLine, OPT_EFFECT, "0.5"));
		}

//		if(cmdLine.hasOption(OPT_POLY_EFFECT_SORT_LONG))
//		{
//			simuMPArgs.setPolyEffectSort();
//		}

		if(cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			simuMPArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		simuMPArgs.setFreq(parseDoubleOptionValue(cmdLine, OPT_FREQ, "0.5"));

		if(cmdLine.hasOption(OPT_UNIF_FREQ_LONG))
		{
			simuMPArgs.setUnifFreq();
		}

		if(cmdLine.hasOption(OPT_FREQ_RANGE_LONG))
		{
			simuMPArgs.setFreqRange(cmdLine.getOptionValues(OPT_FREQ_RANGE_LONG));
		}

		if(cmdLine.hasOption(OPT_FREQ_FILE_LONG))
		{
			simuMPArgs.setFreqFile(cmdLine.getOptionValue(OPT_FREQ_FILE_LONG));
		}

		simuMPArgs.setLD(parseDoubleOptionValue(cmdLine, OPT_LD, "0.0"));

		if(cmdLine.hasOption(OPT_RAND_LD_LONG))
		{
			simuMPArgs.setRandLD();
		}

		if(cmdLine.hasOption(OPT_LD_RANGE_LONG))
		{
			simuMPArgs.setLDRange(cmdLine.getOptionValues(OPT_LD_RANGE_LONG));
		}

		if(cmdLine.hasOption(OPT_HSQ))
		{
			simuMPArgs.setHsq(cmdLine.getOptionValues(OPT_HSQ));
			
		}
		
		if (cmdLine.hasOption(OPT_CM))
		{
			simuMPArgs.setCM(cmdLine.getOptionValue(OPT_CM));
		}

		if (cmdLine.hasOption(OPT_CME))
		{
			simuMPArgs.setCME(cmdLine.getOptionValue(OPT_CME));
		}

		if(cmdLine.hasOption(OPT_MAKE_BED_LONG))
		{
			simuMPArgs.setMakeBed();
		}

		simuMPArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, "2014"));
		return simuMPArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		return new SimulationMPCommandImpl();
	}

	private static final String OPT_SAMPLE_SIZE = "n";
	private static final String OPT_SAMPLE_SIZE_LONG = "sample-size";
	private static final String OPT_SAMPLE_SIZE_LONG_DESC = "Sample size";

	private static final String OPT_MARKER = "m";
	private static final String OPT_MARKER_LONG = "marker";
	private static final String OPT_MARKER_LONG_DESC = "Number of markers.";

	//effect
	private static final String OPT_EFFECT = "effect";
	private static final String OPT_EFFECT_DESC = "Equal effects, 1 by default.";

	private static final String OPT_POLY_EFFECT_LONG = "poly-effect";
	private static final String OPT_POLY_EFFECT_LONG_DESC = "Polygenic (normal distribution) model.";

//	private static final String OPT_POLY_EFFECT_SORT_LONG = "poly-effect-sort";
//	private static final String OPT_POLY_EFFECT_SORT_LONG_DESC = "Sorted polygenic (normal distribution) model.";

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

	private static final String OPT_COR_MAT_LONG = "cor-mat";
	private static final String OPT_CM = "cm";
	private static final String OPT_COR_MAT_LONG_DESC = "correlation matrix";

	private static final String OPT_COR_MAT_E_LONG = "cor-mat-e";
	private static final String OPT_CME = "cme";
	private static final String OPT_COR_MAT_E_LONG_DESC = "correlation matrix for environment";
	
	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_LONG_DESC = "make-bed";

	private static final String OPT_REP = "rep";
	private static final String OPT_REP_DESC = "Replication for simulation. New error will be added to the same breeding values.";
}
