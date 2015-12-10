package gear.subcommands.qpca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class QPCACommand extends Command
{

	@Override
	public String getName()
	{
		return "mwscore";	
	}

	@Override
	public String getDescription()
	{
		return "Generate principal component";
	}

	@Override
	public void prepareOptions(Options options)
	{
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		return null;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		// TODO Auto-generated method stub
		return null;
	}

}
