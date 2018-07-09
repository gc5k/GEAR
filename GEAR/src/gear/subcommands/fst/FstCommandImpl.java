package gear.subcommands.fst;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FSTCommandImpl extends CommandImpl {

	private FSTCommandArguments fstArgs;
	@Override
	public void execute(CommandArguments cmdArgs) {
		fstArgs = (FSTCommandArguments) cmdArgs;
		System.out.println(fstArgs.getGroupFile());
	}

}