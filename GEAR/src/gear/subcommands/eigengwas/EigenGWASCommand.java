package gear.subcommands.eigengwas;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class EigenGWASCommand extends Command {
	public EigenGWASCommand() {
	}

	public String getName() {
		return "egwas";
	}

	public String getDescription() {
		return "EigenGWAS";
	}

	@SuppressWarnings("static-access")
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArgs().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_GUI_DESC).withLongOpt(OPT_GUI_LONG).create());
		
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

	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		EigenGWASCommandArguments eigenArgs = new EigenGWASCommandArguments();
		parseFileArguments((CommandArguments) eigenArgs, cmdLine);

		parseSampleFilterArguments((CommandArguments) eigenArgs, cmdLine);
		parseFamilyFilterArguments((CommandArguments) eigenArgs, cmdLine);

		parseSNPFilterFileArguments((CommandArguments) eigenArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) eigenArgs, cmdLine);

		parseMAFArguments((CommandArguments) eigenArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) eigenArgs, cmdLine);
		parseGENOArguments((CommandArguments) eigenArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) eigenArgs, cmdLine);
		parseMAFRangeArguments((CommandArguments) eigenArgs, cmdLine);

		parsePhenoFileArguments((CommandArguments) eigenArgs, cmdLine);
		parsePhenoIndexArguments((CommandArguments) eigenArgs, cmdLine);
		
		if (cmdLine.hasOption(OPT_GUI_LONG)) {
			eigenArgs.setGUI();
		}
		return eigenArgs;
	}

	protected CommandImpl createCommandImpl() {
		return new EigenGWASCommandImpl();
	}
	
	private final static String OPT_GUI_LONG = "gui";
	private final static String OPT_GUI_DESC = "GUI";

}