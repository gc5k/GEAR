package gear.subcommands.oath.oathbus;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class OATHBusCommand extends Command {

	@Override
	public String getName() {
		return "oath-bus";
	}

	@Override
	public String getDescription() {
		return "Exhaustive OATH";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options) {
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg()
				.isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).hasArg().isRequired().create(OPT_PHE));
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).hasArgs().isRequired().create(OPT_MPHE));
		options.addOption(OptionBuilder.withDescription(OPT_COVAR_DESC).hasArg().isRequired().create(OPT_COVAR));
		options.addOption(
				OptionBuilder.withDescription(OPT_COVAR_NUMBER_DESC).withLongOpt(OPT_COVAR_NUMBER).hasArgs().create());


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

		options.addOption(
				OptionBuilder.withDescription(OPT_REMOVE_OATH_LONG_DESC).withLongOpt(OPT_REMOVE_OATH_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		OATHBusCommandArguments obArgs = new OATHBusCommandArguments();

		parseFileArguments(obArgs, cmdLine);

		parseMAFArguments((CommandArguments) obArgs, cmdLine);
		parseMAXMAFArguments((CommandArguments) obArgs, cmdLine);
		parseGENOArguments((CommandArguments) obArgs, cmdLine);
		parseZeroVarArguments((CommandArguments) obArgs, cmdLine);

		parsePhenoFileArguments((CommandArguments) obArgs, cmdLine);
		parsePhenoIndexArguments((CommandArguments) obArgs, cmdLine);

		
		obArgs.setCovFile(cmdLine.getOptionValue(OPT_COVAR));

		if (cmdLine.hasOption(OPT_COVAR_NUMBER)) {
			obArgs.setCovNumber(cmdLine.getOptionValues(OPT_COVAR_NUMBER));
		}

		if (cmdLine.hasOption(OPT_REMOVE_OATH_LONG)) {
			obArgs.setRemoveInter(cmdLine.getOptionValue(OPT_REMOVE_OATH_LONG));
		}
		return obArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() {
		// TODO Auto-generated method stub
		return new OATHBusCommandImpl();
	}

	private final static String OPT_COVAR = "covar";
	private final static String OPT_COVAR_DESC = "Specify the covariate file";

	private final static String OPT_COVAR_NUMBER = "covar-number";
	private final static String OPT_COVAR_NUMBER_DESC = "Specify the indices for covariate file";
	
	private final static String OPT_REMOVE_OATH_LONG = "remove-oath";
	private final static String OPT_REMOVE_OATH_LONG_DESC = "Remove intermediate results";
}
