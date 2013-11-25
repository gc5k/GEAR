package gear.simulation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class SimuFamilyCommand extends Command
{

	@Override
	public String getName()
	{
		return "simufam";
	}

	@Override
	public String getDescription()
	{
		return "Simulation of discordant nuclear families";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NUM_FAMS_DESC).withLongOpt(OPT_NUM_FAMS_LONG).hasArg().isRequired().create(OPT_NUM_FAMS));
		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKERS_DESC).withLongOpt(OPT_NUM_MARKERS_LONG).hasArg().isRequired().create(OPT_NUM_MARKERS));
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create(OPT_SEED));
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_DESC).withLongOpt(OPT_MAKE_BED_LONG).create(OPT_MAKE_BED));
		options.addOption(OptionBuilder.withDescription(OPT_OUT_DESC).withLongOpt(OPT_OUT_LONG).hasArg().create(OPT_OUT));
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SimuFamilyCommandArguments cmdArgs = new SimuFamilyCommandArguments();
		parseNumberOfFamilies(cmdArgs, cmdLine);
		parseNumberOfMarkers(cmdArgs, cmdLine);
		parseSeed(cmdArgs, cmdLine);
		cmdArgs.setMakeBed(cmdLine.hasOption(OPT_MAKE_BED));
		cmdArgs.setOutRoot(cmdLine.getOptionValue(OPT_OUT, OPT_OUT_DEFAULT));
		return cmdArgs;
	}
	
	private void parseNumberOfFamilies(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		int numFams = 0;
		
		try
		{
			numFams = Integer.parseInt(cmdLine.getOptionValue(OPT_NUM_FAMS));
		}
		catch (NumberFormatException e)
		{
		}
		
		if (numFams <= 0)
		{
			String msg = "";
			msg += "The value of --" + OPT_NUM_FAMS_LONG + "is invalid: '";
			msg += cmdLine.getOptionValue(OPT_NUM_FAMS) + "' is not a valid positive integer.";
			throw new CommandArgumentException(msg);
		}
		
		cmdArgs.setNumberOfFamilies(numFams);
	}
	
	private void parseNumberOfMarkers(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		int numMarkers = 0;
		
		try
		{
			numMarkers = Integer.parseInt(cmdLine.getOptionValue(OPT_NUM_MARKERS));
		}
		catch (NumberFormatException e)
		{
		}
		
		if (numMarkers <= 0)
		{
			String msg = "";
			msg += "The value of --" + OPT_NUM_MARKERS_LONG + "is invalid: '";
			msg += cmdLine.getOptionValue(OPT_NUM_MARKERS) + "' is not a valid positive integer.";
			throw new CommandArgumentException(msg);
		}
		
		cmdArgs.setNumberOfMarkers(numMarkers);
	}
	
	private void parseSeed(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		if (cmdLine.hasOption(OPT_SEED))
		{
			try
			{
				cmdArgs.setSeed(Long.parseLong(cmdLine.getOptionValue(OPT_SEED)));
			}
			catch (NumberFormatException e)
			{
				String msg = "";
				msg += "The value of --" + OPT_SEED_LONG + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_SEED) + "' is not a valid integer.";
				throw new CommandArgumentException(msg);
			}
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimuFamilyCommandImpl();
	}
	
	private static final char OPT_NUM_FAMS = 'f';
	private static final String OPT_NUM_FAMS_LONG = "num-fams";
	private static final String OPT_NUM_FAMS_DESC = "Specify the number of families";
	
	private static final char OPT_NUM_MARKERS = 'm';
	private static final String OPT_NUM_MARKERS_LONG = "num-markers";
	private static final String OPT_NUM_MARKERS_DESC = "Specify the number of markers";
	
	private static final char OPT_SEED = 's';
	private static final String OPT_SEED_LONG = "seed";
	private static final String OPT_SEED_DESC = "Specify the seed of random number generator";
	
	private static final char OPT_MAKE_BED = 'b';
	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_DESC = "Make .bed, .bim and .fam files";

}
