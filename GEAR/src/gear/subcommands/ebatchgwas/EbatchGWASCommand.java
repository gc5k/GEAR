package gear.subcommands.ebatchgwas;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class EbatchGWASCommand extends Command {

	@Override
	public String getName() {
		return "eigengwas";
	}

	@Override
	public String getDescription() {
		return "all in one for eigengwas";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());

		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).withLongOpt(OPT_EV_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_DOM_DESC).create(OPT_DOM));
		options.addOption(OptionBuilder.withDescription(OPT_EPI_DESC).create(OPT_EPI));
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
		EbatchGWASCommandArguments EbatArgs = new EbatchGWASCommandArguments();

		parseFileArguments((CommandArguments) EbatArgs, cmdLine);

		parseSampleFilterArguments((CommandArguments) EbatArgs, cmdLine);
		parseFamilyFilterArguments((CommandArguments) EbatArgs, cmdLine);

		parseSNPFilterFileArguments((CommandArguments) EbatArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) EbatArgs, cmdLine);

		parseMAFArguments((CommandArguments) EbatArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) EbatArgs, cmdLine);
		parseGENOArguments((CommandArguments) EbatArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) EbatArgs, cmdLine);
		parseMAFRangeArguments((CommandArguments) EbatArgs, cmdLine);

		EbatArgs.setEV(parseStringOptionValue(cmdLine, OPT_EV_LONG, "1"));

		if (cmdLine.hasOption(OPT_DOM)) {
			EbatArgs.setDom();
		}
		if (cmdLine.hasOption(OPT_EPI)) {
			EbatArgs.setEpi();
		}

		if (cmdLine.hasOption(OPT_INBRED_LONG)) {
			EbatArgs.setInbred();
		}
		return EbatArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new EbatchGWASCommandImpl();
	}

	private final static String OPT_EV_LONG = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector number";
	private final static String OPT_DOM = "dom";
	private final static String OPT_DOM_DESC = "Specify the dominance effects";
	private final static String OPT_EPI = "epi";
	private final static String OPT_EPI_DESC = "Specify the epistasis effects";

	private final static String OPT_INBRED_LONG = "inbred";
	private final static String OPT_INBRED_DESC = "for inbred line populations.";
}
