package gear.subcommands.oath.nss;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class NSSCommand extends Command {

	@Override
	public String getName() {
		return "nss";
	}

	@Override
	public String getDescription() {
		return "Generates naive summary statistics.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());
//		options.addOption(
//				OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArgs().isRequired().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_DESC).hasArg().isRequired().create(OPT_COVAR));
		options.addOption(
				OptionBuilder.withDescription(OPT_COVAR_NUMBER_DESC).withLongOpt(OPT_COVAR_NUMBER).hasArgs().isRequired().create());

		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

		options.addOption(
				OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());

		options.addOption(OptionBuilder.withDescription(OPT_MAF_DESC).withLongOpt(OPT_MAF_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_MAX_MAF_DESC).withLongOpt(OPT_MAX_MAF_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_GENO_DESC).withLongOpt(OPT_GENO_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_ZERO_VAR_DESC).withLongOpt(OPT_ZERO_VAR_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		NSSCommandArguments nssArgs = new NSSCommandArguments();
		parseFileArguments(nssArgs, cmdLine);
		parseSampleFilterArguments((CommandArguments) nssArgs, cmdLine);
		parseSNPFilterFileArguments((CommandArguments) nssArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) nssArgs, cmdLine);

		parseMAFArguments((CommandArguments) nssArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) nssArgs, cmdLine);
		parseGENOArguments((CommandArguments) nssArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) nssArgs, cmdLine);

		parsePhenoFileArguments((CommandArguments) nssArgs, cmdLine);
		parsePhenoIndexArguments((CommandArguments) nssArgs, cmdLine);

		nssArgs.setCovFile(cmdLine.getOptionValue(OPT_COVAR));
		if (cmdLine.hasOption(OPT_COVAR_NUMBER)) {
			nssArgs.setCovarIndex(cmdLine.getOptionValues(OPT_COVAR_NUMBER));
		}

		return nssArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new NSSCommandImpl();
	}

	protected final static String OPT_COVAR = "covar";
	protected final static String OPT_COVAR_DESC = "Specify the covariate file";

	protected final static String OPT_COVAR_NUMBER = "covar-number";
	protected final static String OPT_COVAR_NUMBER_DESC = "Specify the indices for covariate file";

}
