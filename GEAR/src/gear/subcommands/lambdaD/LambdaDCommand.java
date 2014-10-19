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
		options.addOption(OptionBuilder.withDescription(OPT_META_DESC).withLongOpt(OPT_META_LONG).hasArgs().create(OPT_META));
		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg().create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_DESC).withLongOpt(OPT_META_GZ_LONG).hasArgs().create(OPT_META_GZ));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_BATCH_DESC).withLongOpt(OPT_META_GZ_BATCH_LONG).hasArg().create(OPT_META_GZ_BATCH));

		options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).withLongOpt(OPT_ME_LONG).hasArgs(1).create(OPT_ME));
		options.addOption(OptionBuilder.withDescription(OPT_CC_DESC).hasArgs().create(OPT_CC));
		options.addOption(OptionBuilder.withDescription(OPT_CC_BATCH_LONG_DESC).withLongOpt(OPT_CC_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_QT_DESC).hasArgs().create(OPT_QT));
		options.addOption(OptionBuilder.withDescription(OPT_QT_BATCH_LONG_DESC).withLongOpt(OPT_QT_BATCH_LONG).hasArg().create());
		
		options.addOption(OptionBuilder.withDescription(OPT_KEY_DESC).hasArgs(5).create(OPT_KEY));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_DESC).withLongOpt(OPT_VERBOSE_LONG).create(OPT_VERBOSE));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_GZ_DESC).withLongOpt(OPT_VERBOSE_GZ_LONG).create(OPT_VERBOSE_GZ));
		options.addOption(OptionBuilder.withDescription(OPT_Q_RANGE_DESC).withLongOpt(OPT_Q_RANGE_LONG).hasArgs(2).create(OPT_Q_RANGE));

		options.addOption(OptionBuilder.withDescription(OPT_CM_DESC).hasArg().create(OPT_CM));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		LambdaDCommandArguments lamD = new LambdaDCommandArguments();

		lamD.setMe(OPT_ME_DEFAULT);

		//manual text meta files
		if (cmdLine.hasOption(OPT_META))
		{
			lamD.setMetaFile(cmdLine.getOptionValues(OPT_META));
			lamD.setGZ(false);
			if (cmdLine.hasOption(OPT_QT))
			{
				lamD.setQT(cmdLine.getOptionValues(OPT_QT));
			}
			if (cmdLine.hasOption(OPT_CC))
			{
				lamD.setCC(cmdLine.getOptionValues(OPT_CC));
			}
		}

		//batch text meta batch
		if (cmdLine.hasOption(OPT_META_BATCH))
		{
			lamD.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));
			lamD.setGZ(false);
			if (cmdLine.hasOption(OPT_QT_BATCH_LONG))
			{
				lamD.setQTbatch(cmdLine.getOptionValue(OPT_QT_BATCH_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_BATCH_LONG))
			{
				lamD.setCCbatch(cmdLine.getOptionValue(OPT_CC_BATCH_LONG));
			}
		}

		//me
		if(cmdLine.hasOption(OPT_ME))
		{
			lamD.setMe(cmdLine.getOptionValue(OPT_ME));
		}

		//manual gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ))
		{
			lamD.setMetaFile(cmdLine.getOptionValues(OPT_META_GZ));
			lamD.setGZ(true);
			if (cmdLine.hasOption(OPT_QT))
			{
				lamD.setQT(cmdLine.getOptionValues(OPT_QT));
			}
			if (cmdLine.hasOption(OPT_CC))
			{
				lamD.setCC(cmdLine.getOptionValues(OPT_CC));
			}
		}

		//batch gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ_BATCH))
		{
			lamD.setMetaBatch(cmdLine.getOptionValue(OPT_META_GZ_BATCH));
			lamD.setGZ(true);
			if (cmdLine.hasOption(OPT_QT_BATCH_LONG))
			{
				lamD.setQTbatch(cmdLine.getOptionValue(OPT_QT_BATCH_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_BATCH_LONG))
			{
				lamD.setCCbatch(cmdLine.getOptionValue(OPT_CC_BATCH_LONG));
			}
		}

		if (cmdLine.hasOption(OPT_KEY))
		{
			lamD.setKey(cmdLine.getOptionValues(OPT_KEY));
		}

		if (cmdLine.hasOption(OPT_VERBOSE))
		{
			lamD.setVerbose();
		}
		
		if (cmdLine.hasOption(OPT_VERBOSE_GZ))
		{
			lamD.setVerboseGZ();
		}

		if (cmdLine.hasOption(OPT_Q_RANGE))
		{
			lamD.setQRange(cmdLine.getOptionValues(OPT_Q_RANGE));
		}
		
		if (cmdLine.hasOption(OPT_CM))
		{
			lamD.setCM(cmdLine.getOptionValue(OPT_CM));
		}
		return lamD;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new LambdaDCommandImpl();
	}

	private final static String OPT_META = "m";
	private final static String OPT_META_LONG = "meta";
	private final static String OPT_META_DESC = "The summary statistic files";

	private final static String OPT_META_GZ = "mg";
	private final static String OPT_META_GZ_LONG = "meta-gz";
	private final static String OPT_META_GZ_DESC = "The summary statistic files in gz format";

	private final static String OPT_CC = "cc";
	private final static String OPT_CC_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";

	private final static String OPT_QT = "qt";
	private final static String OPT_QT_DESC = "Quantitative trait: #sample size 1, #sample size 2";


	private final static String OPT_META_BATCH = "mb";
	private final static String OPT_META_BATCH_LONG = "meta-batch";
	private final static String OPT_META_BATCH_DESC = "The summary statistic batch";

	private final static String OPT_META_GZ_BATCH = "mgb";
	private final static String OPT_META_GZ_BATCH_LONG = "meta-gz-batch";
	private final static String OPT_META_GZ_BATCH_DESC = "The summary statistic files in gz format";

	private final static String OPT_CC_BATCH_LONG = "cc-batch";
	private final static String OPT_CC_BATCH_LONG_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";

	private final static String OPT_QT_BATCH_LONG = "qt-batch";
	private final static String OPT_QT_BATCH_LONG_DESC = "Quantitative trait: #sample size 1, #sample size 2";

	private final static String OPT_KEY = "key";
	private final static String OPT_KEY_DESC = "Self defined key workds: snp, beta, se, a1, a2";

	private final static String OPT_VERBOSE = "v";
	private final static String OPT_VERBOSE_LONG = "verbose";
	private final static String OPT_VERBOSE_DESC = "Print test statistic for every pair of meta files.";

	private final static String OPT_VERBOSE_GZ = "vg";
	private final static String OPT_VERBOSE_GZ_LONG = "verbose-gz";
	private final static String OPT_VERBOSE_GZ_DESC = "Print test statistic in gz format for every pair of meta files.";

	private final static String OPT_Q_RANGE = "q";
	private final static String OPT_Q_RANGE_LONG = "q-range";
	private final static String OPT_Q_RANGE_DESC = "cut off values for filtering snps.";

	private final static String OPT_ME = "m";
	private final static String OPT_ME_LONG = "me";
	private final static String OPT_ME_DEFAULT = "30000";
	private final static String OPT_ME_DESC = "effective number of markers.";
	
	private final static String OPT_CM = "cm";
	private final static String OPT_CM_DESC = "correlation matrix";

}
