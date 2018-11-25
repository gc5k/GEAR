package gear.subcommands.grmstat;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class GRMStatCommand extends Command {

	@Override
	public String getName() {
		return "gstat";
	}

	@Override
	public String getDescription() {
		return "Ne and Me";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(
				OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_GRM_TEXT_DESC).withLongOpt(OPT_GRM_TEXT_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_GRM_GZ_DESC).withLongOpt(OPT_GRM_GZ_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		GRMStatCommandArguments gsArgs = new GRMStatCommandArguments();
		parseGRMArguments(gsArgs, cmdLine);
		return gsArgs;
	}

	private void parseGRMArguments(GRMStatCommandArguments gsArgs, CommandLine cmdLine) throws CommandArgumentException {
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TEXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM_GZ_LONG);

		int numFiles = 0;

		if (grmBin != null) {
			gsArgs.setGrmBin(grmBin + ".grm.bin");
			gsArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}

		if (grmText != null) {
			gsArgs.setGrmText(grmText + ".grm");
			gsArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}

		if (grmGZ != null) {
			gsArgs.setGrmGZ(grmGZ + ".grm.gz");
			gsArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}

		if (numFiles == 0) {
			throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", " + OPT_GRM_TEXT_LONG + " or --" + OPT_GRM_GZ_LONG + " must be set.");
		}

		if (numFiles > 1) {
			throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TEXT_LONG + " and --" + OPT_GRM_GZ_LONG + " can be set.");
		}			

	}

	@Override
	protected CommandImpl createCommandImpl() {
		// TODO Auto-generated method stub
		return new GRMStatCommandImpl();
	}

	private final static String OPT_GRM_BIN_LONG = "grm-bin";
	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";

	private final static String OPT_GRM_TEXT_LONG = "grm-txt";
	private final static String OPT_GRM_TEXT_DESC = "Specify the .grm and .grm.id files";

	private final static String OPT_GRM_GZ_LONG = "grm";
	private final static String OPT_GRM_GZ_DESC = "Specify the .grm.gz and .grm.id files";
}
