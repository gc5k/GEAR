package gear.subcommands.metapc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class MetaPCCommand extends Command {

	public MetaPCCommand() {
		addAlias("mpc");
	}

	@Override
	public String getName() {
		return "metapc";
	}

	@Override
	public String getDescription() {
		return "Generating principal components using allele frequency";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg()
				.create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_BATCH_DESC).withLongOpt(OPT_META_GZ_BATCH_LONG)
				.hasArg().create(OPT_META_GZ_BATCH));

		// options.addOption(OptionBuilder.withDescription(OPT_ME_DESC).hasArgs(1).create(OPT_ME));
		// options.addOption(OptionBuilder.withDescription(OPT_NO_WEIGHT_LONG_DESC).withLongOpt(OPT_NO_WEIGHT_LONG).create(OPT_NO_WEIGHT));

		// options.addOption(OptionBuilder.withDescription(OPT_ME_FRAC_LONG_DESC).withLongOpt(OPT_ME_FRAC_LONG).hasOptionalArg().create());

		options.addOption(
				OptionBuilder.withDescription(OPT_CC_SIZE_LONG_DESC).withLongOpt(OPT_CC_SIZE_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_QT_SIZE_LONG_DESC).withLongOpt(OPT_QT_SIZE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_KEY_DESC).hasArgs().create(OPT_KEY));
		options.addOption(OptionBuilder.withDescription(OPT_BETA_DESC).withLongOpt(OPT_BETA_LONG).create(OPT_BETA));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_ATGC_DESC).withLongOpt(OPT_KEEP_ATGC_LONG)
				.create(OPT_KEEP_ATGC));

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		MetaPCCommandArguments metaPC = new MetaPCCommandArguments();

		// metaPC.setMe(OPT_ME_DEFAULT);

		// batch text meta batch
		if (cmdLine.hasOption(OPT_META_BATCH)) {
			metaPC.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));
			metaPC.setGZ(false);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG)) {
				metaPC.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG)) {
				metaPC.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
			}
		}

		if (cmdLine.hasOption(OPT_BETA)) {
			metaPC.setBeta();
		}
		// me
		// if(cmdLine.hasOption(OPT_ME))
		// {
		// metaPC.setMe(cmdLine.getOptionValue(OPT_ME));
		// }
		//
		// if(cmdLine.hasOption(OPT_NO_WEIGHT))
		// {
		// metaPC.setNoWeight();
		// }
		// batch gzip meta files
		if (cmdLine.hasOption(OPT_META_GZ_BATCH)) {
			metaPC.setMetaBatch(cmdLine.getOptionValue(OPT_META_GZ_BATCH));
			metaPC.setGZ(true);
			if (cmdLine.hasOption(OPT_QT_SIZE_LONG)) {
				metaPC.setQTbatch(cmdLine.getOptionValue(OPT_QT_SIZE_LONG));
			}
			if (cmdLine.hasOption(OPT_CC_SIZE_LONG)) {
				metaPC.setCCbatch(cmdLine.getOptionValue(OPT_CC_SIZE_LONG));
			}
		}

		if (cmdLine.hasOption(OPT_KEY)) {
			metaPC.setKey(cmdLine.getOptionValues(OPT_KEY));
		}

		if (cmdLine.hasOption(OPT_KEEP_ATGC)) {
			metaPC.keepATGC();
		}
		// if (cmdLine.hasOption(OPT_VERBOSE))
		// {
		// metaPC.setVerbose();
		// }

		if (cmdLine.hasOption(OPT_CHR_LONG)) {
			metaPC.setChr(cmdLine.getOptionValues(OPT_CHR_LONG));
		}

		// metaPC.setMeFrac(parseDoubleOptionValueInRange(cmdLine, OPT_ME_FRAC_LONG,
		// OPT_ME_FRAC_LONG_DEFAULT, 0.01, 1));

		return metaPC;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new MetaPCCommandImpl();
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

	private final static String OPT_BETA = "b";
	private final static String OPT_BETA_LONG = "beta";
	private final static String OPT_BETA_DESC = "Analysis on genetic effects.";

	// private final static String OPT_VERBOSE = "v";
	// private final static String OPT_VERBOSE_LONG = "verbose";
	// private final static String OPT_VERBOSE_DESC = "Print test statistic for
	// every pair of meta files.";

	// private final static String OPT_ME = "me";
	// private final static String OPT_ME_DEFAULT = "30000";
	// private final static String OPT_ME_DESC = "effective number of markers.";
	//
	// private final static String OPT_NO_WEIGHT = "nw";
	// private final static String OPT_NO_WEIGHT_LONG = "no-weight";
	// private final static String OPT_NO_WEIGHT_LONG_DESC = "Set the sample size
	// same, say 500, for each cohort when calculating Fst.";
	//
	// private final static String OPT_ME_FRAC_LONG = "me-frac";
	// private final static String OPT_ME_FRAC_LONG_DEFAULT = "0.05";
	// private final static String OPT_ME_FRAC_LONG_DESC = "Fraction of minumal
	// markers required, at least 5% by default.";

	private final static String OPT_KEEP_ATGC = "k";
	private final static String OPT_KEEP_ATGC_LONG = "keep-atgc";
	private final static String OPT_KEEP_ATGC_DESC = "Keep parlindormic loci.";
}
