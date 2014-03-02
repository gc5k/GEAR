package gear.subcommands.metawatchdog.decrypt;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class MetaWatchdog2Command extends Command
{
	public MetaWatchdog2Command()
	{
		addAlias("mw2");
	}
	
	@Override
	public String getName()
	{
		return "meta-watchdog2";
	}

	@Override
	public String getDescription()
	{
		return "Find similar individuals between two data sets";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_DATASET1_DESC).withLongOpt(OPT_DATASET1_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_DATASET2_DESC).withLongOpt(OPT_DATASET2_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_CUTOFF_DESC).withLongOpt(OPT_CUTOFF_LONG).hasArg().create(OPT_CUTOFF));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		MetaWatchdog2CommandArguments cmdArgs = new MetaWatchdog2CommandArguments();
		cmdArgs.setDataset1(cmdLine.getOptionValue(OPT_DATASET1_LONG));
		cmdArgs.setDataset2(cmdLine.getOptionValue(OPT_DATASET2_LONG));
		parseCutoff(cmdArgs, cmdLine);
		return cmdArgs;
	}
	
	private void parseCutoff(MetaWatchdog2CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		boolean throwException = false;
		String cutoffString = cmdLine.getOptionValue(OPT_CUTOFF, OPT_CUTOFF_DEFAULT);
		try
		{
			float cutoff = Float.parseFloat(cutoffString);
			if (cutoff < 0 || cutoff > 1)
			{
				throwException = true;
			}
			else
			{
				cmdArgs.setCutoff(cutoff);
			}
		}
		catch (NumberFormatException e)
		{
			throwException = true;
		}
		
		if (throwException)
		{
			String msg = "";
			msg += "'" + cutoffString + "' is not a valid cutoff value.";
			msg += "A cutoff value must be a floating point number ";
			msg += "which is no smaller than 0 and no larger than 1";
			throw new CommandArgumentException(msg);
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new MetaWatchdog2CommandImpl();
	}

	private final static String OPT_DATASET1_LONG = "set1";
	private final static String OPT_DATASET1_DESC = "The score file of the first data set";
	
	private final static String OPT_DATASET2_LONG = "set2";
	private final static String OPT_DATASET2_DESC = "The score file of the second data set";
	
	private final static char OPT_CUTOFF = 'c';
	private final static String OPT_CUTOFF_LONG = "cutoff";
	private final static String OPT_CUTOFF_DEFAULT = "0.95";
	private final static String OPT_CUTOFF_DESC = "Cutoff value, should be no smaller than 0 and no larger than 1, default to " + OPT_CUTOFF_DEFAULT;
}
