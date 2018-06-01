package gear.subcommands.projectedpc;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ProjectedPCCommandImpl extends CommandImpl {

	private ProjectedPCCommandArguments proArgs;

	@Override
	public void execute(CommandArguments cmdArgs) {
		proArgs = (ProjectedPCCommandArguments) cmdArgs;
		System.out.println("haha");
	}

}
