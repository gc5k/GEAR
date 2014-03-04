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
		return "dogpower";
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
		options.addOption(OptionBuilder.withDescription(OPT_ALPHA_DESC).withLongOpt(OPT_ALPHA_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BETA_DESC).withLongOpt(OPT_BETA_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_TESTS_DESC).withLongOpt(OPT_TESTS_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_H2_DESC).withLongOpt(OPT_H2_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		DogPowerCommandArguments cmdArgs = new DogPowerCommandArguments();
		cmdArgs.setAlpha(parseDoubleOptionValueInRange(cmdLine, OPT_ALPHA_LONG, OPT_ALPHA_DEFAULT, 0, 1));
		cmdArgs.setBeta(parseDoubleOptionValueInRange(cmdLine, OPT_BETA_LONG, OPT_BETA_DEFAULT, 0, 1));
		cmdArgs.setTests((int)parseLongOptionValue(cmdLine, OPT_TESTS_LONG, OPT_TESTS_DEFAULT));
		cmdArgs.setH2(parseDoubleOptionValueInRange(cmdLine, OPT_H2_LONG, OPT_H2_DEFAULT, 0, 1));
		cmdArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT));
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new DogPowerCommandImpl();
	}
	
	private static final String OPT_ALPHA_LONG = "alpha";
	private static final String OPT_ALPHA_DEFAULT = "0.05";
	private static final String OPT_ALPHA_DESC = "Specify type I error rate, default to " + OPT_ALPHA_DEFAULT;
	
	private static final String OPT_BETA_LONG = "beta";
	private static final String OPT_BETA_DEFAULT = "0.95";
	private static final String OPT_BETA_DESC = "Specify type II error rate, default to " + OPT_BETA_DEFAULT;
	
	private static final String OPT_TESTS_LONG = "tests";
	private static final String OPT_TESTS_DEFAULT = "100";
	private static final String OPT_TESTS_DESC = "Specify the number of statistical tests, default to " + OPT_TESTS_DEFAULT;
	
	private static final String OPT_H2_LONG = "h2";
	private static final String OPT_H2_DEFAULT = "0.95";
	private static final String OPT_H2_DESC = "Specify the heritability, default to " + OPT_H2_DEFAULT;
}
