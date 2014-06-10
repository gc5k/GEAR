package gear.subcommands.lambdaD;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LambdaDCommand extends Command
{
	public LambdaDCommand()
	{
		addAlias("lam");
	}

	@Override
	public String getName()
	{
		return "lambdad";
	}

	@Override
	public String getDescription()
	{
		return "Calculate Lambda deflation parameter for a pair of meta-analysis";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_META1_DESC).withLongOpt(OPT_META1_LONG).hasArg().create(OPT_META1));
		options.addOption(OptionBuilder.withDescription(OPT_META2_DESC).withLongOpt(OPT_META2_LONG).hasArg().create(OPT_META2));
		options.addOption(OptionBuilder.withDescription(OPT_CC_DESC).hasArgs(4).create(OPT_CC));
		options.addOption(OptionBuilder.withDescription(OPT_QT_DESC).hasArgs(2).create(OPT_QT));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		LambdaDCommandArguments lamD = new LambdaDCommandArguments();

		lamD.setMeta1(cmdLine.getOptionValue(OPT_META1));
		lamD.setMeta2(cmdLine.getOptionValue(OPT_META2));

		if (cmdLine.hasOption(OPT_QT))
		{
			lamD.setQT(cmdLine.getOptionValues(OPT_QT));
		}

		if (cmdLine.hasOption(OPT_CC))
		{
			lamD.setCC(cmdLine.getOptionValues(OPT_CC));
		}
		return lamD;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new LambdaDCommandImpl();
	}

	private final static String OPT_META1 = "m1";
	private final static String OPT_META1_LONG = "meta1";
	private final static String OPT_META1_DESC = "The summary statistics for the first data";
	private final static String OPT_META2 = "m2";
	private final static String OPT_META2_LONG = "meta2";
	private final static String OPT_META2_DESC = "The summary statistics for the second data";

	private final static String OPT_CC = "cc";
	private final static String OPT_CC_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";

	private final static String OPT_QT = "qt";
	private final static String OPT_QT_DESC = "Quantitative trait: #sample size 1, #sample size 2";

}
