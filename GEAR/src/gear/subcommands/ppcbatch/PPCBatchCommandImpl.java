package gear.subcommands.ppcbatch;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.exsnp.ExSNPCommandArguments;
import gear.subcommands.exsnp.ExSNPCommandImpl;
import gear.subcommands.profile.ProfileCommandArguments;
import gear.subcommands.profile.ProfileCommandImpl;
import gear.util.Logger;

public class PPCBatchCommandImpl extends CommandImpl 
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		PPCBatchCommandArguments pbCmdArgs = (PPCBatchCommandArguments) cmdArgs;
		if (!pbCmdArgs.isGreedy())
		{
			ExSNPCommandArguments exSNPcmdArgs = pbCmdArgs.getExSNPCommandArguments();
			exSNPcmdArgs.setOutRoot(pbCmdArgs.getOutRoot());
			ExSNPCommandImpl exSNPImpl = new ExSNPCommandImpl();
			exSNPImpl.execute(exSNPcmdArgs);
			Logger.printUserLog("");
		}

		for(int i = 0; i < pbCmdArgs.getBedFiles().size(); i++)
		{
			ProfileCommandArguments profCommandArguments = pbCmdArgs.getProfileCommandArguments();
			profCommandArguments.setBFile(pbCmdArgs.getBedFiles().get(i));
			if(!pbCmdArgs.isGreedy())
			{
				String scoreExtract = new String(pbCmdArgs.getOutRoot() + ".comsnp");
				profCommandArguments.setIsExtract(scoreExtract);
			}
			profCommandArguments.setResultFile(pbCmdArgs.getOutRoot() + "." + pbCmdArgs.getBedFiles().get(i));
			ProfileCommandImpl profImpl = new ProfileCommandImpl();
			profImpl.execute(profCommandArguments);
			Logger.printUserLog("");

		}
	}
}
