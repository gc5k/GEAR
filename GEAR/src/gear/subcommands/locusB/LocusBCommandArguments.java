package gear.subcommands.locusB;

import gear.subcommands.CommandArguments;

public class LocusBCommandArguments extends CommandArguments {

	public void setInbred() {
		isInbred = true;
	}
	
	public boolean isInbred() {
		return isInbred;
	}

	private boolean isInbred = false;

}
