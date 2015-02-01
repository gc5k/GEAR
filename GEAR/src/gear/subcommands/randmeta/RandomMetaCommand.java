package gear.subcommands.randmeta;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class RandomMetaCommand extends Command
{
	public RandomMetaCommand()
	{
		addAlias("rand");
	}

	@Override
	public String getName()
	{
		return "rand";
	}

	@Override
	public String getDescription()
	{
		return "Calculate correlation for a pair of meta-analysis under the random model.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_META_DESC).withLongOpt(OPT_META_LONG).hasArgs().create(OPT_META));
		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg().create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_DESC).withLongOpt(OPT_META_GZ_LONG).hasArgs().create(OPT_META_GZ));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_BATCH_DESC).withLongOpt(OPT_META_GZ_BATCH_LONG).hasArg().create(OPT_META_GZ_BATCH));

		options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).hasArgs(1).create(OPT_ME));
		options.addOption(OptionBuilder.withDescription(OPT_CC_DESC).hasArgs().create(OPT_CC));
		options.addOption(OptionBuilder.withDescription(OPT_CC_SIZE_LONG_DESC).withLongOpt(OPT_CC_SIZE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_QT_DESC).hasArgs().create(OPT_QT));
		options.addOption(OptionBuilder.withDescription(OPT_QT_SIZE_LONG_DESC).withLongOpt(OPT_QT_SIZE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_KEY_DESC).hasArgs().create(OPT_KEY));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_DESC).withLongOpt(OPT_VERBOSE_LONG).create(OPT_VERBOSE));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_GZ_DESC).withLongOpt(OPT_VERBOSE_GZ_LONG).create(OPT_VERBOSE_GZ));

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
		options.addOption(OptionBuilder.withDescription(OPT_FRQ_DESC).create(OPT_FRQ));

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		RandomMetaCommandArguments randD = new RandomMetaCommandArguments();

		randD.setMe(OPT_ME_DEFAULT);

		//manual text meta files
		if (cmdLine.hasOption(OPT_META))
		{
			randD.setMetaFile(cmdLine.getOptionValues(OPT_META));
			randD.setGZ(false);
			if (cmdLine.hasOption(OPT_QT))
			{
				randD.setQT(cmdLine.getOptionValues(OPT_QT));
			}
			if (cmdLine.hasOption(OPT_CC))
			{
				randD.setCC(cmdLine.getOptionValues(OPT_CC));
			}
		}

		//batch text meta batch
		if (cmdLine.hasOption(OPT_META_BATCH))
		{
			randD.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));
			randD.setGZ(false);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG))
			{
				randD.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG))
			{
				randD.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
			}
		}

		//me
		if(cmdLine.hasOption(OPT_ME))
		{
			randD.setMe(cmdLine.getOptionValue(OPT_ME));
		}

		//manual gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ))
		{
			randD.setMetaFile(cmdLine.getOptionValues(OPT_META_GZ));
			randD.setGZ(true);
			if (cmdLine.hasOption(OPT_QT))
			{
				randD.setQT(cmdLine.getOptionValues(OPT_QT));
			}
			if (cmdLine.hasOption(OPT_CC))
			{
				randD.setCC(cmdLine.getOptionValues(OPT_CC));
			}
		}

		//batch gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ_BATCH))
		{
			randD.setMetaBatch(cmdLine.getOptionValue(OPT_META_GZ_BATCH));
			randD.setGZ(true);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG))
			{
				randD.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG))
			{
				randD.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
			}
		}

		if (cmdLine.hasOption(OPT_KEY))
		{
			randD.setKey(cmdLine.getOptionValues(OPT_KEY));
		}

		if (cmdLine.hasOption(OPT_VERBOSE))
		{
			randD.setVerbose();
		}
		
		if (cmdLine.hasOption(OPT_VERBOSE_GZ))
		{
			randD.setVerboseGZ();
		}

		if (cmdLine.hasOption(OPT_CHR))
		{
			randD.setChr(cmdLine.getOptionValue(OPT_CHR));
		}
		
		if (cmdLine.hasOption(OPT_FRQ))
		{
			randD.setFrq();
		}

		return randD;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new RandomMetaCommandImpl();
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

	private final static String OPT_CC_SIZE_LONG = "cc-size";
	private final static String OPT_CC_SIZE_LONG_DESC = "Case-control study: #case 1, #ctrl 1, #case 2, #ctrl 2";

	private final static String OPT_QT_SIZE_LONG = "qt-size";
	private final static String OPT_QT_SIZE_LONG_DESC = "Quantitative trait: #sample size 1, #sample size 2";

	private final static String OPT_KEY = "key";
	private final static String OPT_KEY_DESC = "Self defined key workds: snp, beta, se, a1, a2, chr, bp, p";

	private final static String OPT_VERBOSE = "v";
	private final static String OPT_VERBOSE_LONG = "verbose";
	private final static String OPT_VERBOSE_DESC = "Print test statistic for every pair of meta files.";

	private final static String OPT_VERBOSE_GZ = "vg";
	private final static String OPT_VERBOSE_GZ_LONG = "verbose-gz";
	private final static String OPT_VERBOSE_GZ_DESC = "Print test statistic in gz format for every pair of meta files.";

	private final static String OPT_ME = "me";
	private final static String OPT_ME_DEFAULT = "30000";
	private final static String OPT_ME_DESC = "effective number of markers.";
	
	private final static String OPT_CHR = "chr";
	private final static String OPT_CHR_DESC = "Choose chromosome for analysis";

	private final static String OPT_FRQ = "frq";
	private final static String OPT_FRQ_DESC = "Using frequency";

}
