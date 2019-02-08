package gear.subcommands.testbed;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.family.plink.PLINKBinaryParser;

public final class TestBedCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs) {
		PLINKBinaryParser parser = new PLINKBinaryParser(cmdArgs.getBed(), cmdArgs.getBim(), cmdArgs.getFam());
		parser.parse();
	}

}
