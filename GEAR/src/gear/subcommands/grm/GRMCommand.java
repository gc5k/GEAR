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

	@Override
	public void prepareOptions(Options options)
	{
	    options.addOption(OptionBuilder.withDescription("Specify PLINK format .bed, .bim and .fam files").withLongOpt("bfile").hasArg().isRequired().create());
	    options.addOption(OptionBuilder.withDescription("Specify the chromosomes for analysis").hasArg().create(OPT_CHR));
	    options.addOption(OptionBuilder.withDescription("Make gz format").create(OPT_GZ));
	    options.addOption(OptionBuilder.withDescription("Make txt format").create(OPT_TXT));
	    options.addOption(OptionBuilder.withDescription("Adjustment for variance").withLongOpt(OPT_VAR_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		GRMArguments grmArgs = new GRMArguments();
	    parseFileArguments(grmArgs, cmdLine);
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
	    	grmArgs.setVar();
	    }
		return grmArgs;
	}

	private void parseFileArguments(GRMArguments grmArgs, CommandLine cmdLine) throws CommandArgumentException 
	{
		String bfile = cmdLine.getOptionValue("bfile");
		String file = cmdLine.getOptionValue("file");

		if ((bfile == null) && (file == null))
		{
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		if ((bfile != null) && (file != null))
		{
			throw new CommandArgumentException("--bfile and --file cannot be set together.");
		}

		grmArgs.setBFile(bfile);
		grmArgs.setFile(file);

		if (cmdLine.hasOption("chr"))
		{
			grmArgs.setChr(cmdLine.getOptionValue("chr"));
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new GRMImpl();
	}

	private static final String OPT_CHR = "chr";
	private static final String OPT_GZ = "gz";
	private static final String OPT_TXT = "txt";
	private static final String OPT_VAR_LONG = "adj-var";

}
