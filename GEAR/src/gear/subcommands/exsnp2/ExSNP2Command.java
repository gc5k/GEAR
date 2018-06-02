package gear.subcommands.exsnp2;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ExSNP2Command extends Command {
	public ExSNP2Command() {
		addAlias("exsnp");
	}

	@Override
	public String getName() {
		return "exsnp";
	}

	@Override
	public String getDescription() {
		return "Extract common snps between files";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_LONG).hasArg().isRequired().create(OPT_BFILE_LONG));
		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().isRequired().create(OPT_BATCH));
//		options.addOption(OptionBuilder.withDescription(OPT_BFILES_DESC).hasArgs().create(OPT_BFILES));
		
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

		ExSNP2CommandArguments esArgs = new ExSNP2CommandArguments();
		parseFileArguments((CommandArguments) esArgs, cmdLine);
		esArgs.setBatch(cmdLine.getOptionValue(OPT_BATCH));

		parseSampleFilterArguments((CommandArguments) esArgs, cmdLine);
		parseFamilyFilterArguments((CommandArguments) esArgs, cmdLine);

		parseSNPFilterFileArguments((CommandArguments) esArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) esArgs, cmdLine);

		parseMAFArguments((CommandArguments) esArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) esArgs, cmdLine);
		parseGENOArguments((CommandArguments) esArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) esArgs, cmdLine);

		parseMAFRangeArguments((CommandArguments) esArgs, cmdLine);

//		if (cmdLine.hasOption(OPT_BFILES)) {
//			esArgs.setBFiles(cmdLine.getOptionValues(OPT_BFILES));
//		}
		return esArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new ExSNP2CommandImpl();
	}

	private String OPT_BATCH = "batch";
	private String OPT_BATCH_DESC = "the batch file for bfiles.";
}
