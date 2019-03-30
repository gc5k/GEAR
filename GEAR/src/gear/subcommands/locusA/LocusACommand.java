package gear.subcommands.locusA;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LocusACommand extends Command {

	@Override
	public String getName() {
		return "locusA";
	}

	@Override
	public String getDescription() {
		return "Statistics for loci";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {

		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());

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

		options.addOption(
				OptionBuilder.withDescription(OPT_INBRED_DESC).withLongOpt(OPT_INBRED_LONG).create());
		
		options.addOption(
				OptionBuilder.withDescription(OPT_THREAD_NUM_LONG_DESC).withLongOpt(OPT_THREAD_NUM_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_THREAD_GREEDY_LONG_DESC).withLongOpt(OPT_THREAD_GREEDY_LONG).hasArg().create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		LocusACommandArguments locArgs = new LocusACommandArguments();

		parseFileArguments(locArgs, cmdLine);

		parseSampleFilterArguments((CommandArguments) locArgs, cmdLine);
		parseFamilyFilterArguments((CommandArguments) locArgs, cmdLine);
		
		parseSNPFilterFileArguments((CommandArguments) locArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) locArgs, cmdLine);

		parseMAFArguments((CommandArguments) locArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) locArgs, cmdLine);
		parseGENOArguments((CommandArguments) locArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) locArgs, cmdLine);
		
		parseMAFRangeArguments((CommandArguments) locArgs, cmdLine);
		parseThreadNumArguments((CommandArguments) locArgs, cmdLine);
		
		if (cmdLine.hasOption(OPT_INBRED_LONG)) {
			locArgs.setInbred();
		}
		return locArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new LocusACommandImpl();
	}

	private static final String OPT_INBRED_LONG = "inbred";
	private static final String OPT_INBRED_DESC = "inbred lines.";
}
