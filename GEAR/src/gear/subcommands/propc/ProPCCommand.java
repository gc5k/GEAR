package gear.subcommands.propc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.profile.ProfileCommand;
import gear.subcommands.profile.ProfileCommandArguments;

public class ProPCCommand extends Command
{

	public ProPCCommand()
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

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().isRequired().create(OPT_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_GREEDY_DESC).create(OPT_GREEDY));
		profCommand.prepareOptions(options);
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		ProPCCommandArguments pcArgs = new ProPCCommandArguments();
		pcArgs.setBatch(cmdLine.getOptionValue(OPT_BATCH));
		pcArgs.setGreedy(cmdLine.hasOption(OPT_GREEDY));
		System.out.println("here1");
		pcArgs.setProfileCommandArguments((ProfileCommandArguments)profCommand.parse(cmdLine));
		System.out.println("here2");
		return pcArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ProPCCommandImpl();
	}
	
	private final static String OPT_BATCH = "batch";
	private final static String OPT_BATCH_DESC = "The batch for all genotype files.";

	private final static String OPT_GREEDY = "greedy";
	private final static String OPT_GREEDY_DESC = "Using all markers.";

	private ProfileCommand profCommand = new ProfileCommand();
}
