package gear.subcommands.grm;

import gear.subcommands.CommandArguments;

public class GRMCommandArguments extends CommandArguments {
	private boolean isGZ = true;
	private boolean isVar = false;
	private boolean isDom = false;
	private boolean isInbred = false;

	public void setGZ() {
		isGZ = true;
	}

	public void setTxt() {
		isGZ = false;
	}

	public boolean isGZ() {
		return isGZ;
	}

	public void setAdjVar() {
		isVar = true;
	}

	public boolean isAdjVar() {
		return isVar;
	}

	public void setDom() {
		isDom = true;
	}

	public boolean isDom() {
		return isDom;
	}

	public void setInbred() {
		isInbred = true;
	}

	public boolean isInbred() {
		return isInbred;
	}
}
