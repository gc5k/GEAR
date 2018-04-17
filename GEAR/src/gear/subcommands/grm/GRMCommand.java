package gear.subcommands.grm;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class GRMCommand extends Command
{

	@Override
	public String getName()
	{
		return "grm";
	}

	@Override
	public String getDescription()
	{
		return "Constructing whole-genome ibs matrix";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
	    options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription("Make gz format").create(OPT_GZ));
	    options.addOption(OptionBuilder.withDescription("Make txt format").create(OPT_TXT));
	    options.addOption(OptionBuilder.withDescription("Dominance").create(OPT_DOM));
	    options.addOption(OptionBuilder.withDescription("MAF").hasArg().create(OPT_MAF));

	    options.addOption(OptionBuilder.withDescription("Adjustment for variance").withLongOpt(OPT_VAR_LONG).create());

	    options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	    options.addOption(OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		GRMArguments grmArgs = new GRMArguments();
	    parseFileArguments((CommandArguments) grmArgs, cmdLine);
	    parseSampleFilterArguments((CommandArguments) grmArgs, cmdLine);
	    parseSNPFilterArguments((CommandArguments) grmArgs, cmdLine);

	    if(cmdLine.hasOption(OPT_GZ))
	    {
	    		grmArgs.setGZ();
	    }
	    if(cmdLine.hasOption(OPT_TXT))
	    {
	    		grmArgs.setTxt();
	    }
	    if(cmdLine.hasOption(OPT_VAR_LONG))
	    {
	    		grmArgs.setAdjVar();
	    }
		if (cmdLine.hasOption(OPT_DOM))
		{
			grmArgs.setDom();
		}
		if (cmdLine.hasOption(OPT_MAF))
		{
			grmArgs.setMAF(cmdLine.getOptionValue(OPT_MAF));
		}
		return grmArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new GRMImpl();
	}

	private static final String OPT_GZ = "gz";
	private static final String OPT_TXT = "txt";
	private static final String OPT_VAR_LONG = "adj-var";
	private static final String OPT_DOM = "dom";
	private static final String OPT_MAF = "maf";

}
