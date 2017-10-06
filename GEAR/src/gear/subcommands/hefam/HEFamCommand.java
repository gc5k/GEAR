package gear.subcommands.hefam;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;

public class HEFamCommand extends Command
{

	@Override
	public String getName()
	{
		return "hefam";
	}

	@Override
	public String getDescription()
	{
		return "Haseman-Elston regression for family-based study";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SD_DESC).create(OPT_SD));
		options.addOption(OptionBuilder.withDescription(OPT_SS_DESC).create(OPT_SS));
		options.addOption(OptionBuilder.withDescription(OPT_CP_DESC).create(OPT_CP));

		options.addOption(OptionBuilder.withDescription(OPT_IBD_DESC).hasArg().create(OPT_IBD));
		options.addOption(OptionBuilder.withDescription(OPT_IBD_LIST_LONG_DESC).withLongOpt(OPT_IBD_LIST_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_COVAR_DESC).hasArg().create(OPT_COVAR));
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_NUMBER_DESC).withLongOpt(OPT_COVAR_NUMBER).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArg().create(OPT_MPHE));

		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).hasArg().create(OPT_KEEP));
		options.addOption(OptionBuilder.withDescription(OPT_SCALE_DESC).create(OPT_SCALE));
		options.addOption(OptionBuilder.withDescription(OPT_JACKKNIFE_DESC).withLongOpt(OPT_JACKKNIFE_LONG).create(OPT_JACKKNIFE));
		options.addOption(OptionBuilder.withDescription(OPT_IBD_CUTOFF_DESC).withLongOpt(OPT_IBD_CUTOFF_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		HEFamCommandArguments heFamArgs = new HEFamCommandArguments();
		parseIBDArguments(heFamArgs, cmdLine);
		heFamArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));

		if (cmdLine.hasOption(OPT_MPHE))
		{
			heFamArgs.setPhenotypeIdx(cmdLine.getOptionValue(OPT_MPHE));			
		}

		if (cmdLine.hasOption(OPT_SCALE))
		{
			heFamArgs.setScale(true);
		}

		if (cmdLine.hasOption(OPT_COVAR))
		{
			heFamArgs.setCovFile(cmdLine.getOptionValue(OPT_COVAR));
		}
		
		if (cmdLine.hasOption(OPT_COVAR_NUMBER))
		{
			heFamArgs.setCovNumber(cmdLine.getOptionValues(OPT_COVAR_NUMBER));
		}
		
		if (cmdLine.hasOption(OPT_KEEP))
		{
			heFamArgs.setKeepFile(cmdLine.getOptionValue(OPT_KEEP));
		}

		if (cmdLine.hasOption(OPT_SD))
		{
			heFamArgs.setSD();
		}

		if (cmdLine.hasOption(OPT_CP))
		{
			heFamArgs.setCP();
		}

		if (cmdLine.hasOption(OPT_SS))
		{
			heFamArgs.setSS();
		}

		if (cmdLine.hasOption(OPT_IBD_CUTOFF_LONG))
		{
			heFamArgs.setIBDcutoff(cmdLine.getOptionValue(OPT_IBD_CUTOFF_LONG));
		}

		if (cmdLine.hasOption(OPT_JACKKNIFE_LONG))
		{
			heFamArgs.setJackknife();
		}

		return heFamArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new HEFamCommandImpl();
	}

	private void parseIBDArguments(HEFamCommandArguments heFamArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		if (cmdLine.hasOption(OPT_IBD))
		{
			String ibd = cmdLine.getOptionValue(OPT_IBD);
			heFamArgs.setIBDGZ(ibd + ".ibd.gz");
			heFamArgs.setIBDID(ibd + ".ibd.id");
			return;
		}

		if (cmdLine.hasOption(OPT_IBD_LIST_LONG))
		{
			heFamArgs.setIBDList(cmdLine.getOptionValue(OPT_IBD_LIST_LONG));
			return;
		}

		Logger.printUserError("No IBD files are provided. One of --" + OPT_IBD + " or --" + OPT_IBD_LIST_LONG + " must be set.");
		System.exit(0);
	}

	private final static String OPT_SD = "sd";
	private final static String OPT_SD_DESC = "Squared difference for a pair of phenotypes";
	private final static String OPT_CP = "cp";
	private final static String OPT_CP_DESC = "Cross-product for a pair for phenotypes";
	private final static String OPT_SS = "ss";
	private final static String OPT_SS_DESC = "Sum of square for a pair for phenotypes";

	private final static String OPT_IBD = "ibd";
	private final static String OPT_IBD_DESC = "Specify the .ibd.gz and .ibd.id files";

	private final static String OPT_IBD_LIST_LONG = "ibd-list";
	private final static String OPT_IBD_LIST_LONG_DESC = "Specify the .ibd.gz and .ibd.id files in list";

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

	private final static String OPT_IBD_CUTOFF_LONG = "grm-cutoff";
	private final static String OPT_IBD_CUTOFF_DESC = "Specify the grm cutfoo for the analysis";

	private final static String OPT_JACKKNIFE_LONG = "jackknife";
	private final static String OPT_JACKKNIFE = "jk";
	private final static String OPT_JACKKNIFE_DESC = "Jackknife";

}
