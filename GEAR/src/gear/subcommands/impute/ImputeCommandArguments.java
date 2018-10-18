package gear.subcommands.impute;

import gear.subcommands.CommandArguments;

public class ImputeCommandArguments extends CommandArguments {

	public void setInbred() {
		isInbred = true;
	}

	public boolean isInbred() {
		return isInbred;
	}

	private boolean isInbred = false;
}
