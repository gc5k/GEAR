package gear.subcommands.eigengwas;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class EigenGWASCommand extends Command
{
	private static final String OPT_PHE_LONG = "pheno";
	private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
	private static final String OPT_MPHE_LONG = "mpheno";
	private static final String OPT_MPHE_DESC = "Specify the phenotype index";
	private static final String OPT_CHR = "chr";
	private static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";

	public EigenGWASCommand()
	{
		addAlias("egwas");
	}

	public String getName()
	{
		return "eigengwas";
	}

	public String getDescription()
	{
		return "Eigen GWAS";
	}

	@SuppressWarnings("static-access")
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).withLongOpt(OPT_PHE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).withLongOpt(OPT_MPHE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
	}

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EigenGWASArguments eigenArgs = new EigenGWASArguments();
		parseFileArguments(eigenArgs, cmdLine);
		eigenArgs.setPhentypeIndex(parseIntOptionValue(cmdLine, "mpheno", "1"));
		eigenArgs.setPhenotypeFile(cmdLine.getOptionValue("pheno"));
		return eigenArgs;
	}

	private void parseFileArguments(EigenGWASArguments eigenArgs, CommandLine cmdLine) throws CommandArgumentException
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

		eigenArgs.setBFile(bfile);
		eigenArgs.setFile(file);

		if (cmdLine.hasOption(OPT_CHR))
		{
			eigenArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}
	}

	protected CommandImpl createCommandImpl()
	{
		return new EigenGWASImpl();
	}
}