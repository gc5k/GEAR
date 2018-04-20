package gear.subcommands.wgrm;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class WGRMCommand extends Command {

	@Override
	public String getName() {
		return "wgrm";
	}

	@Override
	public String getDescription() {
		return "GRM with weights.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_GZ_DESC).create(OPT_GZ));
		options.addOption(OptionBuilder.withDescription(OPT_TXT_DESC).create(OPT_TXT));
		options.addOption(OptionBuilder.withDescription(OPT_DOM_DESC).create(OPT_DOM));

		options.addOption(OptionBuilder.withDescription(OPT_ADJ_VAR_DESC).withLongOpt(OPT_ADJ_VAR_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_WEIGHT_DESC).hasArg().create(OPT_WEIGHT));
		options.addOption(OptionBuilder.withDescription(OPT_WVAR_DESC).create(OPT_WVAR));

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

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		WGRMCommandArguments wgrmArgs = new WGRMCommandArguments();
		parseFileArguments((CommandArguments) wgrmArgs, cmdLine);
		parseSampleFilterArguments((CommandArguments) wgrmArgs, cmdLine);
		parseSNPFilterFileArguments((CommandArguments) wgrmArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) wgrmArgs, cmdLine);

		parseMAFArguments((CommandArguments) wgrmArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) wgrmArgs, cmdLine);
		parseGENOArguments((CommandArguments) wgrmArgs, cmdLine);

		if (cmdLine.hasOption(OPT_GZ)) {
			wgrmArgs.setGZ();
		}

		if (cmdLine.hasOption(OPT_TXT)) {
			wgrmArgs.setTxt();
		}

		if (cmdLine.hasOption(OPT_ADJ_VAR_LONG)) {
			wgrmArgs.setAdjVar();
		}

		if (cmdLine.hasOption(OPT_CHR_LONG)) {
			wgrmArgs.setChr(cmdLine.getOptionValues(OPT_CHR_LONG));
		}

		if (cmdLine.hasOption(OPT_DOM)) {
			wgrmArgs.setDom();
		}

		if (cmdLine.hasOption(OPT_WVAR)) {
			wgrmArgs.setVanRaden();
		}

		if (cmdLine.hasOption(OPT_WEIGHT)) {
			wgrmArgs.setWeightFile(cmdLine.getOptionValue(OPT_WEIGHT));
		}

		return wgrmArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new WGRMCommandImpl();
	}

	private static final String OPT_GZ = "gz";
	private static final String OPT_GZ_DESC = "Make gz format";

	private static final String OPT_TXT = "txt";
	private static final String OPT_TXT_DESC = "Make txt format";

	private static final String OPT_ADJ_VAR_LONG = "adj-var";
	private static final String OPT_ADJ_VAR_DESC = "adjust with variance";

	private static final String OPT_DOM = "dom";
	private static final String OPT_DOM_DESC = "dominance matrix";

	private static final String OPT_WVAR = "vandem";
	private static final String OPT_WVAR_DESC = "RanVandem";

	private static final String OPT_WEIGHT = "weight";
	private static final String OPT_WEIGHT_DESC = "Specify weight file";
}
