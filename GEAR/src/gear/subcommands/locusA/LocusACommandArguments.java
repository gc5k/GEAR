package gear.subcommands.locusA;

import gear.subcommands.CommandArguments;

public class LocusACommandArguments extends CommandArguments {

	public void setInbred() {
		isInbred = true;
	}
	
	public boolean isInbred() {
		return isInbred;
	}

	private boolean isInbred = false;

}
