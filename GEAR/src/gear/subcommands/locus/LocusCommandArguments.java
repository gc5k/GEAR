package gear.subcommands.locus;

import gear.subcommands.CommandArguments;

public class LocusCommandArguments extends CommandArguments {

	public void setInbred() {
		isInbred = true;
	}
	
	public boolean isInbred() {
		return isInbred;
	}

	private boolean isInbred = false;

}
