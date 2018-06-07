package gear.subcommands.fst;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FSTCommand extends Command {

	@Override
	public String getName() {
		return "fst";
	}

	@Override
	public String getDescription() {
		return "Estimating fst";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_GROUP_DESC).withLongOpt(OPT_GROUP_LONG).hasArg().isRequired().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		FSTCommandArguments fstArgs = new FSTCommandArguments();

		parseFileArguments((CommandArguments) fstArgs, cmdLine);
		fstArgs.setGroup(cmdLine.getOptionValue(OPT_GROUP_LONG));
		return fstArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new FSTCommandImpl();
	}

	private static String OPT_GROUP_LONG = "group";
	private static String OPT_GROUP_DESC = "group";
	
}
