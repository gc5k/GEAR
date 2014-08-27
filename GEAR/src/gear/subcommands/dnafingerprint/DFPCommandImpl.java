package gear.subcommands.dnafingerprint;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class DFPCommandImpl extends CommandImpl
{
	private DFPCommandArguments dfpComArgs;
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		dfpComArgs = (DFPCommandArguments) cmdArgs;
	}

}
