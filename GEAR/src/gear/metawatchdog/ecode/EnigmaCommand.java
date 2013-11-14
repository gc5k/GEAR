package gear.metawatchdog.ecode;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.Command;
import gear.CommandArgumentException;
import gear.CommandArguments;
import gear.CommandImpl;

public class EnigmaCommand extends Command
{
	@Override
	public String getName()
	{
		return "enigma";
	}

	@Override
	public String getDescription()
	{
		return "Generate a table of random numbers";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_NUM_COLS_DESC).withLongOpt(OPT_NUM_COLS_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MAP_DESC).withLongOpt(OPT_MAP_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_OUT_DESC).withLongOpt(OPT_OUT_LONG).hasArg().create(OPT_OUT));
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EnigmaCommandArguments cmdArgs = new EnigmaCommandArguments();
		parseSeed(cmdArgs, cmdLine);
		parseNumberOfColumns(cmdArgs, cmdLine);
		cmdArgs.setMapFile(cmdLine.getOptionValue(OPT_MAP_LONG));
		cmdArgs.setOutRoot(cmdLine.getOptionValue(OPT_OUT_LONG, OPT_OUT_DEFAULT));
		return cmdArgs;
	}
	
	private void parseSeed(EnigmaCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String sSeed = cmdLine.getOptionValue(OPT_SEED_LONG, OPT_SEED_DEFAULT);
		try
		{
			cmdArgs.setSeed(Long.valueOf(sSeed));
		}
		catch (NumberFormatException e)
		{
			String msg = "";
			msg += "The option value of --" + OPT_SEED_LONG + " is incorrect. ";
			msg += "'" + sSeed + "' is not a valid integer.";
			throw new CommandArgumentException(msg);
		}
	}
	
	private void parseNumberOfColumns(EnigmaCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String sNumCols = cmdLine.getOptionValue(OPT_NUM_COLS_LONG, OPT_NUM_COLS_DEFAULT);
		
		int numCols = 0;
		boolean valid = true;
		try
		{
			numCols = Integer.valueOf(sNumCols);
		}
		catch (NumberFormatException e)
		{
			valid = false;
		}
		
		if (valid && numCols < 0)
		{
			valid = false;
		}
		
		if (valid)
		{
			cmdArgs.setNumberOfColumns(numCols);
		}
		else
		{
			String msg = "";
			msg += "The option value of --" + OPT_NUM_COLS_LONG + " is incorrect. ";
			msg += "'" + sNumCols + "' is not a valid positive integer.";
			throw new CommandArgumentException(msg);
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new EnigmaCommandImpl();
	}
	
	private final static String OPT_SEED_LONG = "ecode";
	private final static String OPT_SEED_DESC = "Random number seed";
	private final static String OPT_SEED_DEFAULT = "2013";
	
	private final static String OPT_NUM_COLS_LONG = "ecol";
	private final static String OPT_NUM_COLS_DESC = "Number of columns";
	private final static String OPT_NUM_COLS_DEFAULT = "5";
	
	private final static String OPT_MAP_LONG = "refallele";
	private final static String OPT_MAP_DESC = "Map file";
}
