package gear.subcommands.oath.pick;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class OathPickCommand extends Command {

	@Override
	public String getName() {
		return "pick";
	}

	@Override
	public String getDescription() {
		return "pick oath";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_PICK_DESC).withLongOpt(OPT_PICK_LONG).hasArg()
				.isRequired().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		OathPickCommandArguments pkArgs = new OathPickCommandArguments();
		pkArgs.setPickList(cmdLine.getOptionValue(OPT_PICK_LONG));
		return pkArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		// TODO Auto-generated method stub
		return new OathPickCommandImpl();
	}

	protected final static String OPT_PICK_LONG = "pick-list";
	protected final static String OPT_PICK_DESC = "pick list for oath";
}
