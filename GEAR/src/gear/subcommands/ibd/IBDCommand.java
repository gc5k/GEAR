package gear.subcommands.ibd;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class IBDCommand extends Command
{

	@Override
	public String getName()
	{
		return "ibd";
	}

	@Override
	public String getDescription()
	{
		return "Calculate IBD for sibpairs";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
	    options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription("Make gz format").create(OPT_GZ));
	    options.addOption(OptionBuilder.withDescription("Make txt format").create(OPT_TXT));

//	    options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
//	    options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	    options.addOption(OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		IBDCommandArguments ibdArgs = new IBDCommandArguments();
	    parseFileArguments((CommandArguments) ibdArgs, cmdLine);
//	    parseSampleFilterArguments((CommandArguments) ibdArgs, cmdLine);
	    parseSNPFilterArguments((CommandArguments) ibdArgs, cmdLine);

	    if (cmdLine.hasOption(OPT_GZ))
	    {
	    		ibdArgs.setGZ();
	    }
	    if (cmdLine.hasOption(OPT_TXT))
	    {
	    		ibdArgs.setTxt();
	    }
		return ibdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new IBDCommandImpl();
	}

	private static final String OPT_GZ = "gz";
	private static final String OPT_TXT = "txt";
}
