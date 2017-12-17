package gear.subcommands.eigengwasepi;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwasepi.EigenGWASEpiCommandArguments;
import gear.subcommands.eigengwasepi.EigenGWASEpiCommandImpl;

public class EigenGWASEpiCommand extends Command
{
	public EigenGWASEpiCommand()
	{
	}

	public String getName()
	{
		return "egwasepi";
	}

	public String getDescription()
	{
		return "Eigen GWAS a+d+epi";
	}

	@SuppressWarnings("static-access")
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArg().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).hasArg().create(OPT_KEEP));
		options.addOption(OptionBuilder.withDescription(OPT_INBRED_DESC).create(OPT_INBRED));
	}

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EigenGWASEpiCommandArguments eigenEpiArgs = new EigenGWASEpiCommandArguments();

		parseFileArguments(eigenEpiArgs, cmdLine);
		eigenEpiArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));
		eigenEpiArgs.setPhentypeIndex(parseIntOptionValue(cmdLine, OPT_MPHE, "1"));
		if (cmdLine.hasOption(OPT_KEEP))
		{
			eigenEpiArgs.setKeepFile(cmdLine.getOptionValue(OPT_KEEP));
		}
		if (cmdLine.hasOption(OPT_CHR))
		{
			eigenEpiArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}

		if (cmdLine.hasOption(OPT_INBRED))
		{
			eigenEpiArgs.setInbred();
		}
		return eigenEpiArgs;
	}

	private void parseFileArguments(EigenGWASEpiCommandArguments eigenArgs, CommandLine cmdLine) throws CommandArgumentException
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
	}

	protected CommandImpl createCommandImpl()
	{
		return new EigenGWASEpiCommandImpl();
	}

	private static final String OPT_INBRED = "inbred";
	private static final String OPT_INBRED_DESC = "Inbred lines";

	private static final String OPT_PHE = "pheno";
	private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
	private static final String OPT_MPHE = "mpheno";
	private static final String OPT_MPHE_DESC = "Specify the phenotype index";
	private static final String OPT_CHR = "chr";
	private static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";
	private static final String OPT_KEEP = "keep";
	private static final String OPT_KEEP_DESC = "Specify the samples for analysis";
}
