package gear.subcommands.exsnp;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ExSNPCommand extends Command
{
	public ExSNPCommand()
	{
		addAlias("cs");
	}

	@Override
	public String getName()
	{
		return "comsnp";
	}

	@Override
	public String getDescription()
	{
		return "Extract common snps between files";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().create(OPT_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_BFILES_DESC).hasArgs().create(OPT_BFILES));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		ExSNPCommandArguments esArgs = new ExSNPCommandArguments();
		
		boolean hasOpt = false;
		if (cmdLine.hasOption(OPT_BATCH))
		{
			esArgs.setBatch(cmdLine.getOptionValue(OPT_BATCH));
			hasOpt = true;
		}
		if (cmdLine.hasOption(OPT_BFILES))
		{
			esArgs.setBFiles(cmdLine.getOptionValues(OPT_BFILES));
			hasOpt = true;
		}
		return esArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ExSNPCommandImpl();
	}

	private String OPT_BFILES = "bfiles";
	private String OPT_BFILES_DESC = "bfiles from which to extract common snps.";
	
	private String OPT_BATCH = "batch";
	private String OPT_BATCH_DESC = "the batch file for bfiles.";
}
