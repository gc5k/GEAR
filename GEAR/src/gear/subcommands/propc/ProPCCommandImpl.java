package gear.subcommands.propc;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ProPCCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		ProPCCommandArguments proPCArgs = (ProPCCommandArguments) cmdArgs;
		ArrayList<String> bfile = proPCArgs.getbFile();
		for(int i=0; i < bfile.size(); i++)
		{
			System.out.println(bfile.get(i) + " ---");
		}
	}

}
