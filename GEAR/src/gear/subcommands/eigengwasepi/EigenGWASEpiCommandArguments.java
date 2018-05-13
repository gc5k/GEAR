package gear.subcommands.eigengwasepi;

import gear.subcommands.CommandArguments;

public class EigenGWASEpiCommandArguments extends CommandArguments 
{
	private boolean inbred = false;

	public void setInbred()
	{
		inbred = true;
	}

	public boolean isInbred()
	{
		return inbred;
	}
}
