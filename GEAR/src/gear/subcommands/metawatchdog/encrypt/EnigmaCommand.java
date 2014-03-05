package gear.subcommands.metawatchdog.encrypt;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.profile.ProfileCommand;
import gear.subcommands.profile.ProfileCommandArguments;

public class EnigmaCommand extends Command
{
	@Override
	public String getName()
	{
		return "mwscore";
	}

	@Override
	public String getDescription()
	{
		return "Generate a table of random numbers";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_ENCODE_DESC).withLongOpt(OPT_ENCODE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_REFALLELE_DESC).withLongOpt(OPT_REFALLELE_LONG).hasArg().isRequired().create());
		profCommand.setIsCalledByEnigma(true);
		profCommand.prepareOptions(options);
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EnigmaCommandArguments cmdArgs = new EnigmaCommandArguments();
		cmdArgs.setEncodeFile(cmdLine.getOptionValue(OPT_ENCODE_LONG));
		cmdArgs.setRefFile(cmdLine.getOptionValue(OPT_REFALLELE_LONG));
		cmdArgs.setProfileCommandArguments((ProfileCommandArguments)profCommand.parse(cmdLine));
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new EnigmaCommandImpl();
	}

	private final static String OPT_ENCODE_LONG = "encode";
	private final static String OPT_ENCODE_DESC = "The .encode file output by dogpower";

	private final static String OPT_REFALLELE_LONG = "refallele";
	private final static String OPT_REFALLELE_DESC = "Map file";

	private ProfileCommand profCommand = new ProfileCommand();
}
