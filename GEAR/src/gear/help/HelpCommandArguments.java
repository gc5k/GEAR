package gear.help;

import org.apache.commons.cli.CommandLine;

import gear.CommandArgumentException;
import gear.CommandArguments;

public final class HelpCommandArguments extends CommandArguments
{
	public HelpCommandArguments(CommandLine cmdLine) throws CommandArgumentException
	{
		super(cmdLine);
		subcmds = new String[cmdLine.getArgs().length];
		System.arraycopy(cmdLine.getArgs(), 0, subcmds, 0, subcmds.length);
	}
	
	public String[] getSubcommands()
	{
		String[] result = new String[subcmds.length];
		System.arraycopy(subcmds, 0, result, 0, subcmds.length);
		return result;
	}
	
	private String[] subcmds;
}
