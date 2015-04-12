package gear.subcommands.sfst;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SFstCommand extends Command
{

	public SFstCommand()
	{
		addAlias("SFst");
	}
	@Override
	public String getName()
	{
		return "sfst";
	}

	@Override
	public String getDescription()
	{
		return "Calculate Fst using summary statistics";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg().create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_BATCH_DESC).withLongOpt(OPT_META_GZ_BATCH_LONG).hasArg().create(OPT_META_GZ_BATCH));

		options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).hasArgs(1).create(OPT_ME));
		options.addOption(OptionBuilder.withDescription(OPT_NO_WEIGHT_LONG_DESC).withLongOpt(OPT_NO_WEIGHT_LONG).create(OPT_NO_WEIGHT));

		options.addOption(OptionBuilder.withDescription(OPT_ME_FRAC_LONG_DESC).withLongOpt(OPT_ME_FRAC_LONG).hasOptionalArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_CC_SIZE_LONG_DESC).withLongOpt(OPT_CC_SIZE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_QT_SIZE_LONG_DESC).withLongOpt(OPT_QT_SIZE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_KEY_DESC).hasArgs().create(OPT_KEY));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_DESC).withLongOpt(OPT_VERBOSE_LONG).create(OPT_VERBOSE));

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
		options.addOption(OptionBuilder.withDescription(OPT_TOP_DESC).hasArg().create(OPT_TOP));
		options.addOption(OptionBuilder.withDescription(OPT_CLEAN_DESC).create(OPT_CLEAN));

		options.addOption(OptionBuilder.withDescription(OPT_TRIM_DESC).hasOptionalArg().create(OPT_TRIM));

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		SFstCommandArguments lamD = new SFstCommandArguments();

		lamD.setMe(OPT_ME_DEFAULT);

		//batch text meta batch
		if (cmdLine.hasOption(OPT_META_BATCH))
		{
			lamD.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));
			lamD.setGZ(false);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG))
			{
				lamD.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG))
			{
				lamD.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
			}
		}

		//me
		if(cmdLine.hasOption(OPT_ME))
		{
			lamD.setMe(cmdLine.getOptionValue(OPT_ME));
		}

		if(cmdLine.hasOption(OPT_NO_WEIGHT))
		{
			lamD.setNoWeight();
		}
		//batch gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ_BATCH))
		{
			lamD.setMetaBatch(cmdLine.getOptionValue(OPT_META_GZ_BATCH));
			lamD.setGZ(true);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG))
			{
				lamD.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG))
			{
				lamD.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
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

		if (cmdLine.hasOption(OPT_CHR))
		{
			lamD.setChr(cmdLine.getOptionValue(OPT_CHR));
		}

		if (cmdLine.hasOption(OPT_TOP))
		{
			lamD.setTop(cmdLine.getOptionValue(OPT_TOP));
		}

		if (cmdLine.hasOption(OPT_CLEAN))
		{
			lamD.setClean();
		}

		if (cmdLine.hasOption(OPT_TRIM))
		{
			lamD.setTrim(parseDoubleOptionValueInRange(cmdLine, OPT_TRIM, OPT_TRIM_DEFAULT, 0, 0.45));
		}

		lamD.setMeFrac(parseDoubleOptionValueInRange(cmdLine, OPT_ME_FRAC_LONG, OPT_ME_FRAC_LONG_DEFAULT, 0.01, 1));

		return lamD;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SFstCommandImpl();
	}

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

	private final static String OPT_ME = "me";
	private final static String OPT_ME_DEFAULT = "30000";
	private final static String OPT_ME_DESC = "effective number of markers.";

	private final static String OPT_NO_WEIGHT = "nw";
	private final static String OPT_NO_WEIGHT_LONG = "no-weight";
	private final static String OPT_NO_WEIGHT_LONG_DESC = "Set the sample size same, say 500, for each cohort when calculating Fst.";
	
	private final static String OPT_ME_FRAC_LONG = "me-frac";
	private final static String OPT_ME_FRAC_LONG_DEFAULT = "0.05";
	private final static String OPT_ME_FRAC_LONG_DESC = "Fraction of minumal markers required, at least 5% by default.";

	private final static String OPT_CHR = "chr";
	private final static String OPT_CHR_DESC = "Choose chromosome for analysis";

	private final static String OPT_TOP = "top";
	private final static String OPT_TOP_DESC = "Top x files as the reference.";

	private final static String OPT_CLEAN = "clean";
	private final static String OPT_CLEAN_DESC = "No detailed result for each pair of cohorts.";

	private final static String OPT_TRIM = "trim";
	private final static String OPT_TRIM_DEFAULT = "0.001";
	private final static String OPT_TRIM_DESC = "trim off the top and the bottom 5% of markers in calculating lambda_meta and its derived statistics.";

}
