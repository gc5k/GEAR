package gear.subcommands.propc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ProPCCommand extends Command
{

	public ProPCCommand ()
	{
		addAlias("ppc");
	}
	@Override
	public String getName()
	{
		return "propc";
	}

	@Override
	public String getDescription()
	{
		return "Generate projected pc";
	}

	@Override
	public void prepareOptions(Options options)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		ProPCCommandArguments pcArgs = new ProPCCommandArguments();
		return null;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ProPCCommandImpl();
	}

}
