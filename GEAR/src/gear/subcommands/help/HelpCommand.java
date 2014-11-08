package gear.subcommands.help;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;


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
	public void prepareOptions(Options options)
	{
	}
	
	@Override
	protected void printOptionsInEffect(CommandLine cmdLine, String subcmd)
	{
	}
	
	@Override
	public CommandArguments parse(CommandLine cmdLine)
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
