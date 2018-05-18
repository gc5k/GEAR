package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;

public class EigenGWASCommandArguments extends CommandArguments 
{
	private boolean isGUI = false;

	public void setGUI() {
		isGUI = true;
	}
	
	public boolean isGUI() {
		return isGUI;
	}
}