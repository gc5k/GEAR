package gear.subcommands.simuqtreal;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SimulationQTRealCommand extends Command
{

	@Override
	public String getName()
	{
		return "simuqt-real";
	}

	@Override
	public String getDescription()
	{
		return "Simulation polygenic model for quantitative traits";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());

		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_DESC).hasArg().create(OPT_EFFECT));
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_POLY_EFFECT_SORT_LONG_DESC).withLongOpt(OPT_POLY_EFFECT_SORT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_REF_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_REF_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));

		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArg().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_LONG_DESC).withLongOpt(OPT_MAKE_BED_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		SimulationQTRealCommandArguments simuQTRealArgs = new SimulationQTRealCommandArguments();
		parseFileArguments(simuQTRealArgs, cmdLine);
		
		if(cmdLine.hasOption(OPT_REP))
		{
			simuQTRealArgs.setRep(cmdLine.getOptionValue(OPT_REP));
		}

		simuQTRealArgs.setPolyEffect();

		if (cmdLine.hasOption(OPT_EFFECT))
		{
			simuQTRealArgs.setPlainEffect(parseDoubleOptionValue(cmdLine, OPT_EFFECT, "0.5"));
		}

		if (cmdLine.hasOption(OPT_POLY_EFFECT_SORT_LONG))
		{
			simuQTRealArgs.setPolyEffectSort();
		}

		if (cmdLine.hasOption(OPT_REF_EFFECT_FILE_LONG))
		{
			simuQTRealArgs.setRefEffectFile(cmdLine.getOptionValue(OPT_REF_EFFECT_FILE_LONG));
		}

		if(cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			simuQTRealArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		simuQTRealArgs.setHsq(parseDoubleOptionValue(cmdLine, OPT_HSQ, "0.5"));

		
		if(cmdLine.hasOption(OPT_MAKE_BED_LONG))
		{
			simuQTRealArgs.setMakeBed();
		}

		simuQTRealArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, "2014"));
		return simuQTRealArgs;
	}

	private void parseFileArguments(SimulationQTRealCommandArguments simuQTRealArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue("bfile");
		String file = cmdLine.getOptionValue("file");

		if ((bfile == null) && (file == null))
		{
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		if ((bfile != null) && (file != null))
		{
			throw new CommandArgumentException("--bfile and --file cannot be set together.");
		}

		simuQTRealArgs.setBFile(bfile);
		simuQTRealArgs.setFile(file);

	}

	
	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimulationQTRealCommandImpl();
	}

	//effect
	private static final String OPT_EFFECT = "effect";
	private static final String OPT_EFFECT_DESC = "Equal effects, 1 by default.";

	private static final String OPT_POLY_EFFECT_LONG = "poly-effect";
	private static final String OPT_POLY_EFFECT_LONG_DESC = "Polygenic (normal distribution) model.";

	private static final String OPT_POLY_EFFECT_SORT_LONG = "poly-effect-sort";
	private static final String OPT_POLY_EFFECT_SORT_LONG_DESC = "Sorted polygenic (normal distribution) model.";

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	private static final String OPT_REF_EFFECT_FILE_LONG = "rafallele-effect-file";
	private static final String OPT_REF_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	//hsq
	private static final String OPT_HSQ = "hsq";
	private static final String OPT_HSQ_DESC = "Heritability.";

	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_LONG_DESC = "make-bed";

	private static final String OPT_REP = "rep";
	private static final String OPT_REP_DESC = "Replication for simulation. New error will be added to the same breeding values.";
}
