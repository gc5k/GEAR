package gear.subcommands.locus;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LocusCommand extends Command {

	@Override
	public String getName() {
		return "locus";
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

		options.addOption(
				OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());

		options.addOption(
				OptionBuilder.withDescription(OPT_INBRED_DESC).withLongOpt(OPT_INBRED_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		LocusArguments locArgs = new LocusArguments();
		parseFileArguments(locArgs, cmdLine);
		parseSampleFilterArguments((CommandArguments) locArgs, cmdLine);
		parseSNPFilterFileArguments((CommandArguments) locArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) locArgs, cmdLine);

		if (cmdLine.hasOption(OPT_INBRED_LONG)) {
			locArgs.setInbred();
		}
		return locArgs;
	}

	private void parseFileArguments(LocusArguments locusArgs, CommandLine cmdLine) throws CommandArgumentException {
		String bfile = cmdLine.getOptionValue("bfile");

		if ((bfile == null)) {
			throw new CommandArgumentException("No genotypes are provided.");
		}

		locusArgs.setBFile(bfile);
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new LocusImpl();
	}

	private static final String OPT_INBRED_LONG = "inbred";
	private static final String OPT_INBRED_DESC = "inbred lines.";
}
