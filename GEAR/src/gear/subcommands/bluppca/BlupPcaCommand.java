package gear.subcommands.bluppca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class BlupPcaCommand extends Command
{
	@Override
	public String getName()
	{
		return "bluppca";
	}

	@Override
	public String getDescription()
	{
		return "Calculate the SNP effects with BLUP-PCA";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_TEXT_DESC).withLongOpt(OPT_GRM_TEXT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_GZ_DESC).withLongOpt(OPT_GRM_GZ_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).withLongOpt(OPT_PHE_LONG).hasArg().isRequired().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		BlupPcaCommandArguments blupArgs = new BlupPcaCommandArguments();
		parseGRMArguments(blupArgs, cmdLine);
		parseFileArguments(blupArgs, cmdLine);
		blupArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE_LONG));
		return blupArgs;
	}
	
	private void parseGRMArguments(BlupPcaCommandArguments blupArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TEXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM_GZ_LONG);
		
		int numFiles = 0;
		
		if (grmBin != null)
		{
			blupArgs.setGrmBin(grmBin + ".grm.bin");
			blupArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}
		
		if (grmText != null)
		{
			blupArgs.setGrmText(grmText + ".grm");
			blupArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}
		
		if (grmGZ != null)
		{
			blupArgs.setGrmGZ(grmGZ + ".grm.gz");
			blupArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}
		
		if (numFiles == 0)
		{
			throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", " + OPT_GRM_TEXT_LONG + " or --" + OPT_GRM_GZ_LONG + " must be set.");
		}
		
		if (numFiles > 1)
		{
			throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TEXT_LONG + " and --" + OPT_GRM_GZ_LONG + " can be set.");
		}
	}
	
	private void parseFileArguments(BlupPcaCommandArguments blupArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue(OPT_BFILE_LONG);
		String file = cmdLine.getOptionValue(OPT_FILE_LONG);
		
		if (bfile == null && file == null)
		{
			throw new CommandArgumentException("No genotypes are provided. Either --" + OPT_BFILE_LONG + " or --" + OPT_FILE_LONG + " must be set.");
		}
		
		if (bfile != null && file != null)
		{
			throw new CommandArgumentException("--" + OPT_BFILE_LONG + " and --" + OPT_FILE_LONG + " cannot be set together.");
		}
		
		blupArgs.setBFile(bfile);
		blupArgs.setFile(file);
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new BlupPcaCommandImpl();
	}
	
	private final static String OPT_GRM_BIN_LONG = "grm-bin";
	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";
	
	private final static String OPT_GRM_TEXT_LONG = "grm";
	private final static String OPT_GRM_TEXT_DESC = "Specify the .grm and .grm.id files";
	
	private final static String OPT_GRM_GZ_LONG = "grm-gz";
	private final static String OPT_GRM_GZ_DESC = "Specify the .grm.gz and .grm.id files";
	
	private final static String OPT_PHE_LONG = "pheno";
	private final static String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
}
