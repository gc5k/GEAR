package gear.subcommands.metawatchdog.encrypt;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class EnigmaCommand extends Command
{
	@Override
	public String getName()
	{
		return "enigma";
	}

	@Override
	public String getDescription()
	{
		return "Generate a table of random numbers";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_ENCODE_DESC).withLongOpt(OPT_ENCODE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_MAP_DESC).withLongOpt(OPT_MAP_LONG).hasArg().isRequired().create());
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EnigmaCommandArguments cmdArgs = new EnigmaCommandArguments();
		cmdArgs.setEncodeFile(cmdLine.getOptionValue(OPT_ENCODE_LONG));
		cmdArgs.setMapFile(cmdLine.getOptionValue(OPT_MAP_LONG));
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new EnigmaCommandImpl();
	}

	private final static String OPT_ENCODE_LONG = "encode";
	private final static String OPT_ENCODE_DESC = "The .encode file output by dogpower";
	
	private final static String OPT_MAP_LONG = "refallele";
	private final static String OPT_MAP_DESC = "Map file";
}
