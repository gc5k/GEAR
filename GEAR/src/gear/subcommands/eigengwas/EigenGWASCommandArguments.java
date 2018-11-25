package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;

public class EigenGWASCommandArguments extends CommandArguments 
{
	private boolean isGUI = false;
	private boolean isTAB = false;

	public void setGUI() {
		isGUI = true;
	}
	
	public boolean isGUI() {
		return isGUI;
	}

	public void setTAB() {
		isTAB = true;
	}
	
	public boolean isTAB() {
		return isTAB;
	}
}
