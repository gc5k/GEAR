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
		PPCBatchCommandArguments pbArgs = (PPCBatchCommandArguments) cmdArgs;
		if (!pbArgs.isGreedy())
		{
			ExSNPCommandArguments exSNPcmdArgs = pbArgs.getExSNPCommandArguments();
			exSNPcmdArgs.setOutRoot(pbArgs.getOutRoot());
			ExSNPCommandImpl exSNPImpl = new ExSNPCommandImpl();
			exSNPImpl.execute(exSNPcmdArgs);
			Logger.printUserLog("");
		}

		for(int i = 0; i < pbArgs.getBedFiles().size(); i++)
		{
			ProfileCommandArguments profArgs = pbArgs.getProfileCommandArguments();
			profArgs.setBFile(pbArgs.getBedFiles().get(i));
			if(!pbArgs.isGreedy())
			{
				String scoreExtract = new String(pbArgs.getOutRoot() + ".comsnp");
				profArgs.setIsExtractScore(scoreExtract);
				if (pbArgs.getKeepFile() != null) {
					profArgs.setKeepFile(pbArgs.getKeepFile());
				} else if (pbArgs.getRemoveFile() != null) {
					profArgs.setRemoveFile(pbArgs.getRemoveFile());
				}
			}
			profArgs.setResultFile(pbArgs.getOutRoot() + "." + pbArgs.getBedFiles().get(i));
			profArgs.setIsWeighted(false);
			ProfileCommandImpl profImpl = new ProfileCommandImpl();
			profImpl.execute(profArgs);
			Logger.printUserLog("");

		}
	}
}
