package gear.subcommands.hpc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class HpcCommand extends Command
{

	public HpcCommand()
	{
		setStopAtNonOption(true);
	}
	
	@Override
	public String getName()
	{
		return "hpc";
	}

	@Override
	public String getDescription()
	{
		return "Generate (and submit) a qsub script";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NAME_DESC).hasArg().withLongOpt(OPT_NAME_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_RAM_DESC).hasArg().withLongOpt(OPT_RAM_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_EMAIL_DESC).hasArg().withLongOpt(OPT_EMAIL_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_SUBMIT_DESC).withLongOpt(OPT_SUBMIT_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		HpcCommandArguments cmdArgs = new HpcCommandArguments();
		
		cmdArgs.setJobName(cmdLine.getOptionValue(OPT_NAME_LONG, OPT_NAME_DEFAULT));
		
		try
		{
			cmdArgs.setRam(Integer.parseInt(cmdLine.getOptionValue(OPT_RAM_LONG, OPT_RAM_DEFAULT)));
		}
		catch (NumberFormatException e)
		{
			throw new CommandArgumentException("--" + OPT_RAM_LONG + " must be an integer");
		}
		
		cmdArgs.setEmail(cmdLine.getOptionValue(OPT_EMAIL_LONG));
		
		cmdArgs.setIsSubmit(cmdLine.hasOption(OPT_SUBMIT_LONG));
		
		cmdArgs.setNonHpcOptionsAndArguments(cmdLine.getArgs());
		
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new HpcCommandImpl();
	}
	
	private static final String OPT_NAME_LONG = "name";
	private static final String OPT_NAME_DEFAULT = "gear";
	private static final String OPT_NAME_DESC = "Job name, default to '" + OPT_NAME_DEFAULT + "'";
	
	private static final String OPT_RAM_LONG = "ram";
	private static final String OPT_RAM_DEFAULT = "10";
	private static final String OPT_RAM_DESC = "Memory required for the job, default to " + OPT_RAM_DEFAULT + " (GB), must be an integer";
	
	private static final String OPT_EMAIL_LONG = "email";
	private static final String OPT_EMAIL_DESC = "Some cluster requires the user to specify the email for a submitted job. This option can specify the email address";
	
	private static final String OPT_SUBMIT_LONG = "submit";
	private static final String OPT_SUBMIT_DESC = "If this option is set, the generated script will be submitted to cluster";

}
