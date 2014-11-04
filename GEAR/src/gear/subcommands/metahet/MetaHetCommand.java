package gear.subcommands.metahet;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class MetaHetCommand extends Command
{
	public MetaHetCommand()
	{
		addAlias("mh");
	}

	@Override
	public String getName()
	{
		return "meta-het";
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
		options.addOption(OptionBuilder.withDescription(OPT_CC_BATCH_LONG_DESC).withLongOpt(OPT_CC_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).hasArg().create(OPT_ME));
		options.addOption(OptionBuilder.withDescription(OPT_CCO_DESC).hasArg().create(OPT_CCO));
		options.addOption(OptionBuilder.withDescription(OPT_QO_DESC).hasArg().create(OPT_QO));
		options.addOption(OptionBuilder.withDescription(OPT_XM_DESC).hasArg().create(OPT_XM));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		MetaHetArguments mhArgs = new MetaHetArguments();

		//me
		mhArgs.setMe(OPT_ME_DEFAULT);
		if (cmdLine.hasOption(OPT_ME))
		{
			mhArgs.setMe(cmdLine.getOptionValue(OPT_ME));
		}

		//xm
		if (cmdLine.hasOption(OPT_XM))
		{
			mhArgs.setXMFile(cmdLine.getOptionValue(OPT_XM));
		}

		if (cmdLine.hasOption(OPT_QO))
		{
			mhArgs.setQMFile(cmdLine.getOptionValue(OPT_QO));
		}

		if (cmdLine.hasOption(OPT_CCO))
		{
			mhArgs.setCCMFile(cmdLine.getOptionValue(OPT_CCO));
		}

		if (cmdLine.hasOption(OPT_CC_BATCH_LONG))
		{
			mhArgs.setCCbatch(cmdLine.getOptionValue(OPT_CC_BATCH_LONG));
		}
		return mhArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new MetaHetImpl();
	}

	private final static String OPT_XM = "xm";
	private final static String OPT_XM_DESC = "x matrix";

	private final static String OPT_ME = "me";
	private final static String OPT_ME_DEFAULT = "30000";
	private final static String OPT_ME_DESC = "effective number of markers";

	private final static String OPT_QO = "qo";
	private final static String OPT_QO_DESC = "Quantitative trait: #overlapping samples";

	private final static String OPT_CCO = "cco";
	private final static String OPT_CCO_DESC = "Case-control study: #overlapping case:controls";

	private final static String OPT_CC_BATCH_LONG = "cc-batch";
	private final static String OPT_CC_BATCH_LONG_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";
}
