package gear.subcommands.projectedpc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ProjectedPCCommand extends Command {

	@Override
	public String getName() {
		return "propc";
	}

	@Override
	public String getDescription() {
		return "projected pc";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().isRequired().create(OPT_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).hasArg().create(OPT_EV));
		options.addOption(OptionBuilder.withDescription(OPT_INBRED_DESC).withLongOpt(OPT_INBRED_LONG).create());
		
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_FAM_DESC).withLongOpt(OPT_KEEP_FAM_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_REMOVE_FAM_DESC).withLongOpt(OPT_REMOVE_FAM_LONG).hasArg().create());

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
		options.addOption(
				OptionBuilder.withDescription(OPT_MAF_RANGE_DESC).withLongOpt(OPT_MAF_RANGE_LONG).hasArgs().create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		ProjectedPCCommandArguments proArgs = new ProjectedPCCommandArguments();
		parseFileArguments((CommandArguments) proArgs, cmdLine);
		proArgs.setBatch(cmdLine.getOptionValue(OPT_BATCH));
		proArgs.setEV(parseIntOptionValue(cmdLine, OPT_EV, "2"));

		if (cmdLine.hasOption(OPT_INBRED_LONG)) {
			proArgs.setInbred();
		}

		parseSampleFilterArguments((CommandArguments) proArgs, cmdLine);
		parseFamilyFilterArguments((CommandArguments) proArgs, cmdLine);

		parseSNPFilterFileArguments((CommandArguments) proArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) proArgs, cmdLine);

		parseMAFArguments((CommandArguments) proArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) proArgs, cmdLine);
		parseGENOArguments((CommandArguments) proArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) proArgs, cmdLine);

		parseMAFRangeArguments((CommandArguments) proArgs, cmdLine);

		return proArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new ProjectedPCCommandImpl();
	}

	private static final String OPT_BATCH = "batch";
	private static final String OPT_BATCH_DESC = "Specify batch file";

	private static final String OPT_EV = "ev";
	private static final String OPT_EV_DESC = "Specify eigenvector number";

	private final static String OPT_INBRED_LONG = "inbred";
	private final static String OPT_INBRED_DESC = "for inbred line populations.";
}
