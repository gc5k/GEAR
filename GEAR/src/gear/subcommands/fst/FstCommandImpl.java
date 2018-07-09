package gear.subcommands.fst;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FstCommandImpl extends CommandImpl {

	private FstCommandArguments fstArgs;
	@Override
	public void execute(CommandArguments cmdArgs) {
		fstArgs = (FstCommandArguments) cmdArgs;
		System.out.println(fstArgs.getGroupFile());
	}

}