package gear.subcommands.impute;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.ibd.IBDCommandArguments;

public class ImputeCommand extends Command {

	@Override
	public String getName() {
		return "impute";
	}

	@Override
	public String getDescription() {
		return "naive imputation";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());

		options.addOption(
				OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
		options.addOption(
				OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		ImputeCommandArguments imputeArgs = new ImputeCommandArguments();
		parseFileArguments((CommandArguments) imputeArgs, cmdLine);
		parseSNPFilterFileArguments((CommandArguments) imputeArgs, cmdLine);
		parseSNPFilterChromosomeArguments((CommandArguments) imputeArgs, cmdLine);
		parseSampleFilterArguments((CommandArguments) imputeArgs, cmdLine);

		return imputeArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		return new ImputeCommandImpl();
	}

}
