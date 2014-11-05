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
		options.addOption(OptionBuilder.withDescription(OPT_CC_SIZE_LONG_DESC).withLongOpt(OPT_CC_SIZE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).hasArg().create(OPT_ME));
		options.addOption(OptionBuilder.withDescription(OPT_OM_DESC).hasArg().create(OPT_OM));
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

		if (cmdLine.hasOption(OPT_OM))
		{
			mhArgs.setOMFile(cmdLine.getOptionValue(OPT_OM));
		}

		if (cmdLine.hasOption(OPT_CC_SIZE_LONG))
		{
			mhArgs.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
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

	private final static String OPT_OM = "om";
	private final static String OPT_OM_DESC = "#overlapping samples. The diagonal are sample size for cohorts, the lower triangle are overlapping samples. For case-control studies, the lower triangle elements are overalpping cases and the upper triangle elements are overlapping controls.";

	private final static String OPT_CC_SIZE_LONG = "cc-size";
	private final static String OPT_CC_SIZE_LONG_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";
}
