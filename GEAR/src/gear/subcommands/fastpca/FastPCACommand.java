package gear.subcommands.fastpca;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.grm.GRMArguments;

public class FastPCACommand extends Command
{
	@Override
	public String getName()
	{
		return "fpca";
	}

	@Override
	public String getDescription()
	{
		return "Calculate FastPCA";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
	    options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());

//		options.addOption(OptionBuilder.withDescription(OPT_GRM_BIN_DESC).withLongOpt(OPT_GRM_BIN_LONG).create());
//		options.addOption(OptionBuilder.withDescription(OPT_GRM_TXT_DESC).withLongOpt(OPT_GRM_TXT_LONG).create());
//		options.addOption(OptionBuilder.withDescription(OPT_GRM_DESC).create(OPT_GRM));
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).hasArg().create(OPT_EV));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).hasArg().create(OPT_KEEP));
		options.addOption(OptionBuilder.withDescription(OPT_PROP_DESC).hasArg().create(OPT_PROP));
	    options.addOption(OptionBuilder.withDescription(OPT_VAR_LONG_DESC).withLongOpt(OPT_VAR_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		FastPCACommandArguments fpcaArgs = new FastPCACommandArguments();
	    parseFileArguments(fpcaArgs, cmdLine);
//	    parseGRMArguments(fpcaArgs, cmdLine);

		if (cmdLine.hasOption(OPT_EV))
		{
			fpcaArgs.setEV(cmdLine.getOptionValue(OPT_EV));
		}
		if (cmdLine.hasOption(OPT_KEEP))
		{
			fpcaArgs.setKeepFile(cmdLine.getOptionValue(OPT_KEEP));
		}
		if (cmdLine.hasOption(OPT_PROP))
		{
			fpcaArgs.setProp(cmdLine.getOptionValue(OPT_PROP));
		}
		if (cmdLine.hasOption(OPT_VAR_LONG))
		{
			fpcaArgs.setVar();
		}
		fpcaArgs.setSeed(parseLongOptionValueInRange(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT, 0, Long.MAX_VALUE));
		return fpcaArgs;
	}

	private void parseFileArguments(FastPCACommandArguments grmArgs, CommandLine cmdLine) throws CommandArgumentException 
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
			grmArgs.setBFile(bfile);			
		}
		if (file != null)
		{
			grmArgs.setFile(file);			
		}

	}

//	private void parseGRMArguments(FastPCACommandArguments fpcaArgs, CommandLine cmdLine) throws CommandArgumentException
//	{
//		String grmBin = cmdLine.getOptionValue(OPT_GRM_BIN_LONG);
//		String grmText = cmdLine.getOptionValue(OPT_GRM_TXT_LONG);
//		String grmGZ = cmdLine.getOptionValue(OPT_GRM);
//		
//		int numFiles = 0;
//		
//		if (grmBin != null)
//		{
//			fpcaArgs.setGrmBin(grmBin + ".grm.bin");
//			fpcaArgs.setGrmID(grmBin + ".grm.id");
//			++numFiles;
//		}
//
//		if (grmText != null)
//		{
//			fpcaArgs.setGrmText(grmText + ".grm");
//			fpcaArgs.setGrmID(grmText + ".grm.id");
//			++numFiles;
//		}
//
//		if (grmGZ != null)
//		{
//			fpcaArgs.setGrmGZ(grmGZ + ".grm.gz");
//			fpcaArgs.setGrmID(grmGZ + ".grm.id");
//			++numFiles;
//		}
//
//		if (numFiles == 0)
//		{
//			throw new CommandArgumentException("No GRM is provided. One of --" + OPT_GRM_BIN_LONG + ", " + OPT_GRM_TXT_LONG + " or --" + OPT_GRM + " must be set.");
//		}
//		
//		if (numFiles > 1)
//		{
//			throw new CommandArgumentException("At most one of --" + OPT_GRM_BIN_LONG + ", --" + OPT_GRM_TXT_LONG + " and --" + OPT_GRM + " can be set.");
//		}
//	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new FastPCACommandImpl();
	}

//	private final static String OPT_GRM_BIN_LONG = "grm-bin";
//	private final static String OPT_GRM_BIN_DESC = "Specify the .grm.bin and .grm.id files";
//
//	private final static String OPT_GRM_TXT_LONG = "grm-txt";
//	private final static String OPT_GRM_TXT_DESC = "Specify the .grm and .grm.id files";
//
//	private final static String OPT_GRM = "grm";
//	private final static String OPT_GRM_DESC = "Specify the .grm.gz and .grm.id files";
	
	private final static String OPT_EV = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector numbers";

	private final static String OPT_KEEP = "keep";
	private final static String OPT_KEEP_DESC = "Specify the samples for analysis";

	private final static String OPT_PROP = "prop";
	private final static String OPT_PROP_DESC = "Specify the proportion";
	
	private static final String OPT_VAR_LONG = "adj-var";
	private static final String OPT_VAR_LONG_DESC = "Adjust the grm with variance.";
}
