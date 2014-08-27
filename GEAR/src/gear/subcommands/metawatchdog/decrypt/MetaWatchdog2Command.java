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
		addAlias("mw");
	}

	@Override
	public String getName()
	{
		return "meta-watchdog";
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
		options.addOption(OptionBuilder.withDescription(OPT_CHISQ_DESC).hasArg().create(OPT_CHISQ));
		options.addOption(OptionBuilder.withDescription(OPT_B_DESC).withLongOpt(OPT_B_LONG).hasArg().create(OPT_B));

		options.addOption(OptionBuilder.withDescription(OPT_ENCODE_DESC).withLongOpt(OPT_ENCODE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_DESC).withLongOpt(OPT_VERBOSE_LONG).create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		MetaWatchdog2CommandArguments cmdArgs = new MetaWatchdog2CommandArguments();
		cmdArgs.setDataset1(cmdLine.getOptionValue(OPT_DATASET1_LONG));
		cmdArgs.setDataset2(cmdLine.getOptionValue(OPT_DATASET2_LONG));
		
		if (cmdLine.hasOption(OPT_CHISQ))
		{
			cmdArgs.setChisq(parseDoubleOptionValueInRange(cmdLine, OPT_CHISQ, OPT_CHISQ_DEFAULT, 0, 10));
		}
		if (cmdLine.hasOption(OPT_B_LONG))
		{
			cmdArgs.setRegB(parseDoubleOptionValueInRange(cmdLine, OPT_B_LONG, OPT_B_DEFAULT, 0, 1));
		}
		if (cmdLine.hasOption(OPT_ENCODE_LONG))
		{
			cmdArgs.setEncodeFile(cmdLine.getOptionValue(OPT_ENCODE_LONG));			
		}
		if (cmdLine.hasOption(OPT_VERBOSE_LONG))
		{
			cmdArgs.setVerbose();
		}

		return cmdArgs;
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

	private final static String OPT_CHISQ = "chisq";
	private final static String OPT_CHISQ_DESC = "chi square test";
	private final static String OPT_CHISQ_DEFAULT = "0.05";

	private final static char OPT_B = 'b';
	private final static String OPT_B_LONG = "reg";
	private final static String OPT_B_DEFAULT = "0.95";
	private final static String OPT_B_DESC = "Cutoff value for regression. value should be no smaller than 0 and no larger than 1, default to " + OPT_B_DEFAULT;
	
	private final static String OPT_ENCODE_LONG = "encode";
	private final static String OPT_ENCODE_DESC = "The .encode file output by dogpower";
	
	private final static String OPT_VERBOSE_LONG = "verbose";
	private final static String OPT_VERBOSE_DESC = "Print out all results.";
}
