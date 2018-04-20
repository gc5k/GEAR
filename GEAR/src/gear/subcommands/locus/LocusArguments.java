package gear.subcommands.locus;

import gear.subcommands.CommandArguments;

public class LocusArguments extends CommandArguments
{

	public void setInbred() {
		isInbred = true;
	}
	
	public boolean isInbred() {
		return isInbred;
	}

	private boolean isInbred = false;
}
