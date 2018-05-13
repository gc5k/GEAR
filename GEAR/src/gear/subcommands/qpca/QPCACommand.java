package gear.subcommands.qpca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class QPCACommand extends Command {

	public QPCACommand() {

	}

	@Override
	public String getName() {
		return "pca";
	}

	@Override
	public String getDescription() {
		return "Calculate PCA";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(
				OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_GRM_TXT_DESC).withLongOpt(OPT_GRM_TXT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_DESC).hasArg().create(OPT_GRM));
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).hasArg().create(OPT_EV));
		// options.addOption(OptionBuilder.withDescription(OPT_VAR_LONG_DESC).withLongOpt(OPT_VAR_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		QPCACommandArguments qpcaArgs = new QPCACommandArguments();
		parseGRMArguments(qpcaArgs, cmdLine);
		qpcaArgs.setEV(parseStringOptionValue(cmdLine, OPT_EV, "10"));
		return qpcaArgs;
	}

	private void parseGRMArguments(QPCACommandArguments qpcaArgs, CommandLine cmdLine) throws CommandArgumentException {
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM);

		int numFiles = 0;

		if (grmBin != null) {
			qpcaArgs.setGrmBin(grmBin + ".grm.bin");
			qpcaArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}

		if (grmText != null) {
			qpcaArgs.setGrmText(grmText + ".grm");
			qpcaArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}

		if (grmGZ != null) {
			qpcaArgs.setGrmGZ(grmGZ + ".grm.gz");
			qpcaArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}

		if (numFiles == 0) {
			throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", "
					+ OPT_GRM_TXT_LONG + " or --" + OPT_GRM + " must be set.");
		}

		if (numFiles > 1) {
			throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TXT_LONG
					+ " and --" + OPT_GRM + " can be set.");
		}
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new QPCACommandImpl();
	}

	private final static String OPT_GRM_BIN_LONG = "grm-bin";
	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";

	private final static String OPT_GRM_TXT_LONG = "grm-txt";
	private final static String OPT_GRM_TXT_DESC = "Specify the .grm and .grm.id files";

	private final static String OPT_GRM = "grm";
	private final static String OPT_GRM_DESC = "Specify the .grm.gz and .grm.id files";

	private final static String OPT_EV = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector numbers";

	//
	// private static final String OPT_VAR_LONG = "adj-var";
	// private static final String OPT_VAR_LONG_DESC = "Adjust the grm with
	// variance.";
}
