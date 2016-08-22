package gear.subcommands.reml;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class REMLCommand extends Command {

	@Override
	public String getName()
	{
		return "reml";
	}

	@Override
	public String getDescription() 
	{
		return "reml";
	}

	@Override
	public void prepareOptions(Options options) 
	{
		options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_TEXT_DESC).withLongOpt(OPT_GRM_TEXT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_GZ_DESC).withLongOpt(OPT_GRM_GZ_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_LIST_DESC).withLongOpt(OPT_GRM_LIST_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_DESC).hasArg().create(OPT_COVAR));
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_NUMBER_DESC).withLongOpt(OPT_COVAR_NUMBER).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArg().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_MINQUE_DESC).create(OPT_MINQUE));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		REMLCommandArguments remlArgs = new REMLCommandArguments();
		parseGRMArguments(remlArgs, cmdLine);
		remlArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));
		if(cmdLine.hasOption(OPT_MPHE))
		{
			remlArgs.setPhenotypeIdx(cmdLine.getOptionValue(OPT_MPHE));			
		}
		if (cmdLine.hasOption(OPT_MINQUE))
		{
			remlArgs.setMINQUE(true);
		}

		if (cmdLine.hasOption(OPT_COVAR))
		{
			remlArgs.setCovFile(cmdLine.getOptionValue(OPT_COVAR));
		}
		
		if (cmdLine.hasOption(OPT_COVAR_NUMBER))
		{
			remlArgs.setCovNumber(cmdLine.getOptionValues(OPT_COVAR_NUMBER));
		}

		return remlArgs;
	}

	private void parseGRMArguments(REMLCommandArguments remlArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TEXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM_GZ_LONG);
		String grmList = cmdLine.getOptionValue(OPT_GRM_LIST_LONG);

		int numFiles = 0;
		
		if (grmBin != null)
		{
			remlArgs.setGrmBin(grmBin + ".grm.bin");
			remlArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}
		
		if (grmText != null)
		{
			remlArgs.setGrmText(grmText + ".grm");
			remlArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}
		
		if (grmGZ != null)
		{
			remlArgs.setGrmGZ(grmGZ + ".grm.gz");
			remlArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}

		if (grmList != null)
		{
			remlArgs.setGrmList(cmdLine.getOptionValue(OPT_GRM_LIST_LONG));
		}
		else
		{
			if (numFiles == 0)
			{
				throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", " + OPT_GRM_TEXT_LONG + " or --" + OPT_GRM_GZ_LONG + " must be set.");
			}
			if (numFiles > 1)
			{
				throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TEXT_LONG + " and --" + OPT_GRM_GZ_LONG + " can be set.");
			}			
		}
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		return new REMLCommandImpl();
	}

	private final static String OPT_GRM_LIST_LONG = "grm-list";
	private final static String OPT_GRM_LIST_DESC = "Specify the .grm files";
	
	private final static String OPT_GRM_BIN_LONG = "grm-bin";
	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";

	private final static String OPT_GRM_TEXT_LONG = "grm-txt";
	private final static String OPT_GRM_TEXT_DESC = "Specify the .grm and .grm.id files";

	private final static String OPT_GRM_GZ_LONG = "grm";
	private final static String OPT_GRM_GZ_DESC = "Specify the .grm.gz and .grm.id files";

	private final static String OPT_PHE = "pheno";
	private final static String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";

	private final static String OPT_MPHE = "mpheno";
	private final static String OPT_MPHE_DESC = "Specify the index for phenotype";

	private final static String OPT_MINQUE = "minque";
	private final static String OPT_MINQUE_DESC = "Minimum norm quadratic unbiased estimation";

	private final static String OPT_COVAR = "covar";
	private final static String OPT_COVAR_DESC = "Specify the covariate file";

	private final static String OPT_COVAR_NUMBER = "covar-number";
	private final static String OPT_COVAR_NUMBER_DESC = "Specify the indices for covariate file";

}
