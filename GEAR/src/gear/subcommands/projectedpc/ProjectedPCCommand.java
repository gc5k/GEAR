package gear.subcommands.projectedpc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ProjectedPCCommand extends Command {

	@Override
	public String getName() {
		return "propc";
	}

	@Override
	public String getDescription() {
		return "projected pc";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().isRequired().create(OPT_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).hasArg().create(OPT_EV));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		ProjectedPCCommandArguments proArgs = new ProjectedPCCommandArguments();
		parseFileArguments((CommandArguments) proArgs, cmdLine);
		proArgs.setEV(parseIntOptionValue(cmdLine, OPT_EV, "2"));

		return proArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new ProjectedPCCommandImpl();
	}

	private static final String OPT_BATCH = "batch";
	private static final String OPT_BATCH_DESC = "Specify batch file";

	private static final String OPT_EV = "ev";
	private static final String OPT_EV_DESC = "Specify eigenvector number";
}
