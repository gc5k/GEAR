package gear.subcommands.lmdr;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LMDRCommand extends Command
{

	@Override
	public String getName()
	{
		return "lmdr";
	}

	@Override
	public String getDescription()
	{
		return "Linear model implementation for MDR";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_CV_DESC).hasArg().create(OPT_CV));
		
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		LMDRCommandArguments lmdrArgs = new LMDRCommandArguments();
		lmdrArgs.setCV(cmdLine.getOptionValue(OPT_CV));
		return lmdrArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new LMDRCommandImpl();
	}

	private String OPT_CV = "cv";
	private String OPT_CV_DESC = "cross-validation";
}
