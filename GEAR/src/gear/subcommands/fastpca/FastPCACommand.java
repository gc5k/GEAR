package gear.subcommands.fastpca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FastPCACommand extends Command {
	@Override
	public String getName() {
		return "fpca";
	}

	@Override
	public String getDescription() {
		return "Calculate FastPCA";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		// options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());

		// options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).create());
		// options.addOption(OptionBuilder.withDescription(OPT_GRM_TXT_DESC).withLongOpt(OPT_GRM_TXT_LONG).create());
		// options.addOption(OptionBuilder.withDescription(OPT_GRM_DESC).create(OPT_GRM));
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).hasArg().create(OPT_EV));

		options.addOption(OptionBuilder.withDescription(OPT_PROP_DESC).hasArg().create(OPT_PROP));
		options.addOption(OptionBuilder.withDescription(OPT_VAR_LONG_DESC).withLongOpt(OPT_VAR_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

		options.addOption(
				OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		FastPCACommandArguments fpcaArgs = new FastPCACommandArguments();
		parseFileArguments((CommandArguments) fpcaArgs, cmdLine);
		parseSampleFilterArguments((CommandArguments) fpcaArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) fpcaArgs, cmdLine);

		if (cmdLine.hasOption(OPT_EV)) {
			fpcaArgs.setEV(cmdLine.getOptionValue(OPT_EV));
		}
		if (cmdLine.hasOption(OPT_PROP)) {
			fpcaArgs.setProp(cmdLine.getOptionValue(OPT_PROP));
		}
		if (cmdLine.hasOption(OPT_VAR_LONG)) {
			fpcaArgs.setAdjVar();
		}
		fpcaArgs.setSeed(parseLongOptionValueInRange(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT, 0, Long.MAX_VALUE));
		return fpcaArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new FastPCACommandImpl();
	}

	private final static String OPT_EV = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector numbers";

	private final static String OPT_PROP = "prop";
	private final static String OPT_PROP_DESC = "Specify the proportion";

	private static final String OPT_VAR_LONG = "adj-var";
	private static final String OPT_VAR_LONG_DESC = "Adjust the grm with variance.";
}
