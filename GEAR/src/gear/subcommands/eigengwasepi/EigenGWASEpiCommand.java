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

	    options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	    options.addOption(OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());

		options.addOption(OptionBuilder.withDescription(OPT_MAF_DESC).withLongOpt(OPT_MAF_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_MAX_MAF_DESC).withLongOpt(OPT_MAX_MAF_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GENO_DESC).withLongOpt(OPT_GENO_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_ZERO_VAR_DESC).withLongOpt(OPT_ZERO_VAR_LONG).create());
	    options.addOption(
				OptionBuilder.withDescription(OPT_MAF_RANGE_DESC).withLongOpt(OPT_MAF_RANGE_LONG).hasArgs().create());

	    options.addOption(OptionBuilder.withDescription(OPT_INBRED_DESC).create(OPT_INBRED));
	}

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EigenGWASEpiCommandArguments eigenEpiArgs = new EigenGWASEpiCommandArguments();

		parseFileArguments(eigenEpiArgs, cmdLine);
		
	    parseSampleFilterArguments((CommandArguments) eigenEpiArgs, cmdLine);
	    parseSNPFilterFileArguments((CommandArguments) eigenEpiArgs, cmdLine);
	    parseSNPFilterChromosomeArguments((CommandArguments) eigenEpiArgs, cmdLine);

		parseMAFArguments((CommandArguments) eigenEpiArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) eigenEpiArgs, cmdLine);
		parseGENOArguments((CommandArguments) eigenEpiArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) eigenEpiArgs, cmdLine);
		parseMAFRangeArguments((CommandArguments) eigenEpiArgs, cmdLine);

		parsePhenoFileArguments((CommandArguments) eigenEpiArgs, cmdLine);
		parsePhenoIndexArguments((CommandArguments) eigenEpiArgs, cmdLine);

		if (cmdLine.hasOption(OPT_INBRED)) eigenEpiArgs.setInbred();
		return eigenEpiArgs;
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
}
