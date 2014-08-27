package gear.subcommands.metawatchdog.powercalculator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class DogPowerCommand extends Command
{
	@Override
	public String getName()
	{
		return "mwpower";
	}

	@Override
	public String getDescription()
	{
		return "Calculate how many score columns are needed for meta-watchdog analysis.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		//generic parameters
		options.addOption(OptionBuilder.withDescription(OPT_ALPHA_DESC).withLongOpt(OPT_ALPHA_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MISSING_DESC).withLongOpt(OPT_MISSING_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_TEST_DESC).withLongOpt(OPT_TEST_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());

		//regression parameters
		options.addOption(OptionBuilder.withDescription(OPT_REGRESSION_DESC).withLongOpt(OPT_REGRESSION_LONG).hasOptionalArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BETA_DESC).withLongOpt(OPT_BETA_LONG).hasArg().create());

		//chisq parameters
		options.addOption(OptionBuilder.withDescription(OPT_CHISQ_DESC).withLongOpt(OPT_CHISQ_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		DogPowerCommandArguments cmdArgs = new DogPowerCommandArguments();
		cmdArgs.setAlpha(parseDoubleOptionValueInRange(cmdLine, OPT_ALPHA_LONG, OPT_ALPHA_DEFAULT, 0, 1));	
		cmdArgs.setMissingRate(parseDoubleOptionValueInRange(cmdLine, OPT_MISSING_LONG, OPT_MISSING_DEFAULT, 0.01, 0.05));

		cmdArgs.setTests(parseLongOptionValue(cmdLine, OPT_TEST_LONG, OPT_TEST_DEFAULT));

		cmdArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT));

		if (cmdLine.hasOption(OPT_REGRESSION_LONG))
		{
			cmdArgs.setRegression(parseDoubleOptionValueInRange(cmdLine, OPT_REGRESSION_LONG, OPT_BVALUE_DEFAULT, 0, 1));
		}
		cmdArgs.setBeta(parseDoubleOptionValueInRange(cmdLine, OPT_BETA_LONG, OPT_BETA_DEFAULT, 0, 1));

		if (cmdLine.hasOption(OPT_CHISQ_LONG))
		{
			cmdArgs.setChisq();
		}
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new DogPowerCommandImpl();
	}

//generic
	private static final String OPT_ALPHA_LONG = "alpha";
	private static final String OPT_ALPHA_DEFAULT = "0.05";
	private static final String OPT_ALPHA_DESC = "Specify type I error rate, default to " + OPT_ALPHA_DEFAULT;

	private static final String OPT_TEST_LONG = "test";
	private static final String OPT_TEST_DEFAULT = "1000";
	private static final String OPT_TEST_DESC = "Specify the number of statistical tests, default to " + OPT_TEST_DEFAULT;

	private static final String OPT_MISSING_LONG = "miss";
	private static final String OPT_MISSING_DEFAULT = "0.015";
	private static final String OPT_MISSING_DESC = "Missing rate used for detecting overlapping individuals, default to " + OPT_MISSING_DEFAULT;

//regression
	private static final String OPT_REGRESSION_LONG = "reg";
	private static final String OPT_BVALUE_DEFAULT = "0.95";
	private static final String OPT_REGRESSION_DESC = "Regression method.";

	private static final String OPT_BETA_LONG = "beta";
	private static final String OPT_BETA_DEFAULT = "0.05";
	private static final String OPT_BETA_DESC = "Specify type II error rate, default to " + OPT_BETA_DEFAULT;

//chisq
	private static final String OPT_CHISQ_LONG = "chisq";
	private static final String OPT_CHISQ_DESC = "Chisq method for detecting overlapping samples.";

}
