package gear.help;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import gear.Command;
import gear.CommandArguments;
import gear.CommandImpl;


public final class HelpCommand extends Command
{
	public HelpCommand()
	{
		addAlias("?");
		addAlias("h");
	}

	@Override
	public String getName()
	{
		return "help";
	}

	@Override
	public String getDescription()
	{
		return "Print help information";
	}

	@Override
	public String getLongDescription()
	{
		return "Describe the usage of this program or its subcommands.";
	}
	
	@Override
	protected void prepareOptions(Options options)
	{
	}
	
	@Override
	protected void printOptionsInEffect(CommandLine cmdLine)
	{
	}
	
	@Override
	protected CommandArguments parse(CommandLine cmdLine)
	{
		HelpCommandArguments cmdArgs = new HelpCommandArguments();
		cmdArgs.setSubcommands(cmdLine.getArgs());
		return cmdArgs;
	}
	
	@Override
	protected CommandImpl createCommandImpl()
	{
		return new HelpCommandImpl();
	}
}
