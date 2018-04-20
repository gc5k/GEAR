package gear.subcommands.eigengwasdom;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class EigenGWASDomCommand extends Command
{
	private static final String OPT_PHE = "pheno";
	private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
	private static final String OPT_MPHE = "mpheno";
	private static final String OPT_MPHE_DESC = "Specify the phenotype index";

	public EigenGWASDomCommand()
	{
	}

	public String getName()
	{
		return "egwasd";
	}

	public String getDescription()
	{
		return "Eigen GWAS a+d";
	}

	@SuppressWarnings("static-access")
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArg().create(OPT_MPHE));

	    options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	    options.addOption(OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());

	}

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EigenGWASDomCommandArguments eigenArgs = new EigenGWASDomCommandArguments();
		parseFileArguments(eigenArgs, cmdLine);
	    parseSampleFilterArguments((CommandArguments) eigenArgs, cmdLine);

	    parseSNPFilterFileArguments((CommandArguments) eigenArgs, cmdLine);
	    parseSNPFilterChromosomeArguments((CommandArguments) eigenArgs, cmdLine);

		eigenArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));
		eigenArgs.setPhentypeIndex(parseIntOptionValue(cmdLine, OPT_MPHE, "1"));

		return eigenArgs;
	}

	protected CommandImpl createCommandImpl()
	{
		return new EigenGWASDomImpl();
	}
}
