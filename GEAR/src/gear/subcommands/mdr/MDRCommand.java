package gear.subcommands.mdr;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class MDRCommand extends Command {
	private static final String OPT_PHE_LONG = "pheno";
	private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
	private static final String OPT_MPHE_LONG = "mpheno";
	private static final String OPT_MPHE_DESC = "Specify the phenotype index";

	private static final String OPT_COV = "cov";
	private static final String OPT_COV_DESC = "Specify the covariates index";

	private static final String OPT_CC = "cc";
	private static final String OPT_CC_DESC = "Case-control design";

	@Override
	public String getName() {
		return "mdr";
	}

	public String getDescription() {
		return "Gene-gene interaction";
	}

	@SuppressWarnings("static-access")
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription("Specify PLINK format .bed, .bim and .fam files")
				.withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_PHE_DESC).withLongOpt(OPT_PHE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).withLongOpt(OPT_MPHE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_COV_DESC).hasArgs().create(OPT_COV));
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_CC_DESC).create(OPT_CC));
	}

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		MDRArguments mdrArgs = new MDRArguments();
		parseFileArguments(mdrArgs, cmdLine);
		mdrArgs.setCovIndexes(cmdLine.getOptionValues(OPT_COV));
		mdrArgs.setPhenotypeIndex(parseIntOptionValue(cmdLine, "mpheno", "1"));
		mdrArgs.setPhenotypeFile(cmdLine.getOptionValue("pheno"));
		mdrArgs.setCC(cmdLine.hasOption(OPT_CC));
		return mdrArgs;
	}

	protected CommandImpl createCommandImpl() {
		return new MDRImpl();
	}
}
