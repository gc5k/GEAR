package gear.subcommands.hereg;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class HERegCommand extends Command 
{

	@Override
	public String getName()
	{
		return "he";
	}

	@Override
	public String getDescription() 
	{
		return "Haseman-Elston regression";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) 
	{
		options.addOption(OptionBuilder.withDescription(OPT_SD_DESC).create(OPT_SD));
		options.addOption(OptionBuilder.withDescription(OPT_SS_DESC).create(OPT_SS));
		options.addOption(OptionBuilder.withDescription(OPT_CP_DESC).create(OPT_CP));

		options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_TEXT_DESC).withLongOpt(OPT_GRM_TEXT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_GZ_DESC).withLongOpt(OPT_GRM_GZ_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_LIST_DESC).withLongOpt(OPT_GRM_LIST_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_DESC).hasArg().create(OPT_COVAR));
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_NUMBER_DESC).withLongOpt(OPT_COVAR_NUMBER).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArg().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).hasArg().create(OPT_KEEP));
		options.addOption(OptionBuilder.withDescription(OPT_SCALE_DESC).create(OPT_SCALE));
		options.addOption(OptionBuilder.withDescription(OPT_JACKKNIFE_DESC).withLongOpt(OPT_JACKKNIFE_LONG).create(OPT_JACKKNIFE));
		options.addOption(OptionBuilder.withDescription(OPT_GRM_CUTOFF_DESC).withLongOpt(OPT_GRM_CUTOFF_LONG).hasArg().create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		HERegCommandArguments heArgs = new HERegCommandArguments();
		parseGRMArguments(heArgs, cmdLine);
		heArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));
		if(cmdLine.hasOption(OPT_MPHE))
		{
			heArgs.setPhenotypeIdx(cmdLine.getOptionValue(OPT_MPHE));			
		}
		if (cmdLine.hasOption(OPT_SCALE))
		{
			heArgs.setScale(true);
		}

		if (cmdLine.hasOption(OPT_COVAR))
		{
			heArgs.setCovFile(cmdLine.getOptionValue(OPT_COVAR));
		}
		
		if (cmdLine.hasOption(OPT_COVAR_NUMBER))
		{
			heArgs.setCovNumber(cmdLine.getOptionValues(OPT_COVAR_NUMBER));
		}
		
		if (cmdLine.hasOption(OPT_KEEP))
		{
			heArgs.setKeepFile(cmdLine.getOptionValue(OPT_KEEP));
		}

		if (cmdLine.hasOption(OPT_SD))
		{
			heArgs.setSD();
		}
		
		if (cmdLine.hasOption(OPT_CP))
		{
			heArgs.setCP();
		}
		
		if (cmdLine.hasOption(OPT_SS))
		{
			heArgs.setSS();
		}
		
		if (cmdLine.hasOption(OPT_GRM_CUTOFF_LONG))
		{
			heArgs.setGRMcutoff(cmdLine.getOptionValue(OPT_GRM_CUTOFF_LONG));
		}
		
		if (cmdLine.hasOption(OPT_JACKKNIFE_LONG))
		{
			heArgs.setJackknife();
		}
		return heArgs;
	}

	private void parseGRMArguments(HERegCommandArguments heArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TEXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM_GZ_LONG);
		String grmList = cmdLine.getOptionValue(OPT_GRM_LIST_LONG);

		int numFiles = 0;
		
		if (grmBin != null)
		{
			heArgs.setGrmBin(grmBin + ".grm.bin");
			heArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}
		
		if (grmText != null)
		{
			heArgs.setGrmText(grmText + ".grm");
			heArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}
		
		if (grmGZ != null)
		{
			heArgs.setGrmGZ(grmGZ + ".grm.gz");
			heArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}

		if (grmList != null)
		{
			heArgs.setGrmList(cmdLine.getOptionValue(OPT_GRM_LIST_LONG));
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
		return new HERegCommandImpl();
	}

	private static String OPT_SD = "sd";
	private static String OPT_SD_DESC = "Squared difference for a pair of phenotypes";
	private static String OPT_CP = "cp";
	private static String OPT_CP_DESC = "Cross-product for a pair for phenotypes";
	private static String OPT_SS = "ss";
	private static String OPT_SS_DESC = "Sum of square for a pair for phenotypes";

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

	private final static String OPT_SCALE = "scale";
	private final static String OPT_SCALE_DESC = "Standardize the phenotype";

	private final static String OPT_COVAR = "covar";
	private final static String OPT_COVAR_DESC = "Specify the covariate file";

	private final static String OPT_COVAR_NUMBER = "covar-number";
	private final static String OPT_COVAR_NUMBER_DESC = "Specify the indices for covariate file";
	
	private final static String OPT_KEEP = "keep";
	private final static String OPT_KEEP_DESC = "Specify the individuals to be used in the analysis";

	private final static String OPT_GRM_CUTOFF_LONG = "grm-cutoff";
	private final static String OPT_GRM_CUTOFF_DESC = "Specify the grm cutfoo for the analysis";
	
	private final static String OPT_JACKKNIFE_LONG = "jackknife";
	private final static String OPT_JACKKNIFE = "jk";
	private final static String OPT_JACKKNIFE_DESC = "Jackknife";

}
