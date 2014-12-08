package gear.subcommands.lmdr;

import gear.subcommands.CommandArguments;

public class LMDRCommandArguments extends CommandArguments
{

	public void setCV(String cv)
	{
		this.cv = Integer.parseInt(cv);
	}
	
	public int getCV()
	{
		return cv;
	}
	private int cv = 5;
}
