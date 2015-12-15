package gear.subcommands.qpca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class QPCACommand extends Command
{
	@Override
	public String getName()
	{
		return "pca";
	}

	@Override
	public String getDescription()
	{
		return "Calculate PCA";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_TEXT_DESC).withLongOpt(OPT_GRM_TEXT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GRM_DESC).withLongOpt(OPT_GRM_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).withLongOpt(OPT_EV_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		QPCACommandArguments qpcaArgs = new QPCACommandArguments();
		parseGRMArguments(qpcaArgs, cmdLine);
		qpcaArgs.setEV(cmdLine.getOptionValue(OPT_EV_LONG));
		return qpcaArgs;
	}

	private void parseGRMArguments(QPCACommandArguments qpcaArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
		String grmText = cmdLine.getOptionValue(OPT_GRM_TEXT_LONG);
		String grmGZ = cmdLine.getOptionValue(OPT_GRM_LONG);
		
		int numFiles = 0;
		
		if (grmBin != null)
		{
			qpcaArgs.setGrmBin(grmBin + ".grm.bin");
			qpcaArgs.setGrmID(grmBin + ".grm.id");
			++numFiles;
		}

		if (grmText != null)
		{
			qpcaArgs.setGrmText(grmText + ".grm");
			qpcaArgs.setGrmID(grmText + ".grm.id");
			++numFiles;
		}

		if (grmGZ != null)
		{
			qpcaArgs.setGrmGZ(grmGZ + ".grm.gz");
			qpcaArgs.setGrmID(grmGZ + ".grm.id");
			++numFiles;
		}

		if (numFiles == 0)
		{
			throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", " + OPT_GRM_TEXT_LONG + " or --" + OPT_GRM_LONG + " must be set.");
		}
		
		if (numFiles > 1)
		{
			throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TEXT_LONG + " and --" + OPT_GRM_LONG + " can be set.");
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new QPCACommandImpl();
	}

	private final static String OPT_GRM_BIN_LONG = "grm-bin";
	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";
	
	private final static String OPT_GRM_TEXT_LONG = "grm-text";
	private final static String OPT_GRM_TEXT_DESC = "Specify the .grm and .grm.id files";
	
	private final static String OPT_GRM_LONG = "grm";
	private final static String OPT_GRM_DESC = "Specify the .grm.gz and .grm.id files";
	
	private final static String OPT_EV_LONG = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector numbers";
}
