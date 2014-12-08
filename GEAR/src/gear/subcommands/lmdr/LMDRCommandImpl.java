package gear.subcommands.lmdr;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LMDRCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		lmdrArgs = (LMDRCommandArguments) cmdArgs;
		System.out.println(lmdrArgs.getCV());
		
	}
	
	private LMDRCommandArguments lmdrArgs;

}
