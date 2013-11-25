package gear.subcommands.help;

import gear.subcommands.CommandArguments;

public final class HelpCommandArguments extends CommandArguments
{
	public void setSubcommands(String[] subcmds)
	{
		this.subcmds = new String[subcmds.length];
		System.arraycopy(subcmds, 0, this.subcmds, 0, subcmds.length);
	}
	
	public String[] getSubcommands()
	{
		String[] result = new String[subcmds.length];
		System.arraycopy(subcmds, 0, result, 0, subcmds.length);
		return result;
	}
	
	private String[] subcmds;
}
