package gear.subcommands.oath.nss;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwas.EigenGWASArguments;

public class NSSCommand extends Command 
{

	@Override
	public String getName() 
	{
		return "nss";
	}

	@Override
	public String getDescription() 
	{
		return "Generates naive summary statistics.";
	}

	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).withLongOpt(OPT_PHE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).withLongOpt(OPT_MPHE_LONG).hasArgs().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		NSSCommandArguments nssArgs = new NSSCommandArguments();
		parseFileArguments(nssArgs, cmdLine);
		nssArgs.setPhentypeIndex(cmdLine.getOptionValues("mpheno"));
		nssArgs.setPhenotypeFile(cmdLine.getOptionValue("pheno"));
		return nssArgs;
	}

	private void parseFileArguments(NSSCommandArguments nssArgs, CommandLine cmdLine) throws CommandArgumentException
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

		nssArgs.setBFile(bfile);
		nssArgs.setFile(file);

		if (cmdLine.hasOption(OPT_CHR))
		{
			nssArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new NSSCommandImpl();
	}

	private static final String OPT_PHE_LONG = "pheno";
	private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
	private static final String OPT_MPHE_LONG = "mpheno";
	private static final String OPT_MPHE_DESC = "Specify the phenotype indicex";
	private static final String OPT_CHR = "chr";
	private static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";

}