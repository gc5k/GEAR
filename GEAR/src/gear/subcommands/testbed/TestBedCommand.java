package gear.subcommands.testbed;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class TestBedCommand extends Command {

	@Override
	public String getName() {
		return "testbed";
	}

	@Override
	public String getDescription() {
		return "Just read bed file to test reading performance";
	}

	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());			
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		TestBedCommandArguments testBedArgs = new TestBedCommandArguments();
		String bfile = cmdLine.getOptionValue(OPT_BFILE_LONG);
		testBedArgs.setBFile(bfile);
		return testBedArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new TestBedCommandImpl();
	}

}
