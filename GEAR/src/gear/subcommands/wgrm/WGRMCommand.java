package gear.subcommands.wgrm;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class WGRMCommand extends Command 
{

	@Override
	public String getName()
	{
		return "wgrm";
	}

	@Override
	public String getDescription()
	{
		return "GRM with weights.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
	    options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
	    options.addOption(OptionBuilder.withDescription(OPT_GZ_DESC).create(OPT_GZ));
	    options.addOption(OptionBuilder.withDescription(OPT_TXT_DESC).create(OPT_TXT));
	    options.addOption(OptionBuilder.withDescription(OPT_DOM_DESC).create(OPT_DOM));
	    options.addOption(OptionBuilder.withDescription(OPT_MAF_DESC).hasArg().create(OPT_MAF));
	    options.addOption(OptionBuilder.withDescription(OPT_VAR_DESC).withLongOpt(OPT_VAR_LONG).create());

	    options.addOption(OptionBuilder.withDescription(OPT_WEIGHT_DESC).hasArg().create(OPT_WEIGHT));
	    options.addOption(OptionBuilder.withDescription(OPT_WVAR_DESC).create(OPT_WVAR));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		WGRMCommandArguments wgrmArgs = new WGRMCommandArguments();
	    parseFileArguments(wgrmArgs, cmdLine);
	    if(cmdLine.hasOption(OPT_GZ))
	    {
	    	wgrmArgs.setGZ();
	    }

	    if(cmdLine.hasOption(OPT_TXT))
	    {
	    	wgrmArgs.setTxt();
	    }

	    if(cmdLine.hasOption(OPT_VAR_LONG))
	    {
	    	wgrmArgs.setAdjVar();
	    }

		if (cmdLine.hasOption(OPT_CHR))
		{
			wgrmArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}

		if (cmdLine.hasOption(OPT_DOM))
		{
			wgrmArgs.setDom();
		}

		if (cmdLine.hasOption(OPT_MAF))
		{
			wgrmArgs.setMAF(cmdLine.getOptionValue(OPT_MAF));
		}

		if (cmdLine.hasOption(OPT_WVAR))
		{
			wgrmArgs.setVanRaden();
		}

		if (cmdLine.hasOption(OPT_WEIGHT))
		{
			wgrmArgs.setWeightFile(cmdLine.getOptionValue(OPT_WEIGHT));
		}

		return wgrmArgs;
	}

	private void parseFileArguments(WGRMCommandArguments wgrmArgs, CommandLine cmdLine) throws CommandArgumentException 
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

		if (bfile != null)
		{
			wgrmArgs.setBFile(bfile);			
		}
		if (file != null)
		{
			wgrmArgs.setFile(file);			
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new WGRMCommandImpl();
	}

	private static final String OPT_CHR = "chr";
	private static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";

	private static final String OPT_GZ = "gz";
	private static final String OPT_GZ_DESC = "Make gz format";

	private static final String OPT_TXT = "txt";
	private static final String OPT_TXT_DESC = "Make txt format";

	private static final String OPT_VAR_LONG = "adj-var";
	private static final String OPT_VAR_DESC = "adjust with variance";

	private static final String OPT_DOM = "dom";
	private static final String OPT_DOM_DESC = "dominance matrix";

	private static final String OPT_MAF = "maf";
	private static final String OPT_MAF_DESC = "Minor allele frequency cutoff";

	private static final String OPT_WVAR = "vandem";
	private static final String OPT_WVAR_DESC = "RanVandem";

	private static final String OPT_WEIGHT= "weight";
	private static final String OPT_WEIGHT_DESC = "Specify weight file";
}
