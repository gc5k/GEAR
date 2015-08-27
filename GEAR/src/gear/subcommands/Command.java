package gear.subcommands;

import gear.AboutInfo;
import gear.util.Logger;

import java.util.Collections;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public abstract class Command implements Comparable<Command>
{
	public abstract String getName();
	
	public boolean hasAlias(String alias)
	{
		return aliases.contains(alias);
	}
	
	protected void addAlias(String alias)
	{
		aliases.add(alias);
	}
	
	public Set<String> getAliases()
	{
		return Collections.unmodifiableSet(aliases);
	}
	
	@Override
	public boolean equals(Object obj)
	{
		String s = null;
		
		if (obj instanceof String)
		{
			s = (String)obj;
		}
		else if (obj instanceof Command)
		{
			s = ((Command)obj).getName();
		}
		
		return s != null && (getName().equals(s) || hasAlias(s));
	}
	
	public int compareTo(Command otherCmd)
	{
		if (equals(otherCmd))
		{
			return 0;
		}
		return getName().compareTo(otherCmd.getName());
	}
	
	private Set<String> aliases = new TreeSet<String>();
	
	public abstract String getDescription();
	
	public String getLongDescription()
	{
		return getDescription();
	}
	
	public String getFullDescription()
	{
		return "";
	}
	
	public abstract void prepareOptions(Options options);
	
	@SuppressWarnings("static-access")
	public Options getOptions()
	{
		Options options = new Options();
		options.addOption(OptionBuilder.withDescription(OPT_OUT_DESC).withLongOpt(OPT_OUT_LONG).hasArg().create(OPT_OUT));
		prepareOptions(options);
		return options;
	}
	
	public abstract CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException;
	protected abstract CommandImpl createCommandImpl();
	
	protected void printOptionsInEffect(CommandLine cmdLine, String subcmd)
	{
		Logger.printUserLog("Subcommand: " + subcmd);
		Logger.printUserLog("Options in effect: ");
		
		@SuppressWarnings("rawtypes")
		Iterator optIter = cmdLine.iterator();
		while (optIter.hasNext())
		{
			Option opt = (Option)optIter.next();
			String line = opt.hasLongOpt() ? "\t--" + opt.getLongOpt() : "\t--" + opt.getOpt();
			String[] argValues = opt.getValues();
			if (argValues != null)
			{
				for (String value : argValues)
				{
					line += " " + value; 
				}
			}
			Logger.printUserLog(line);
		}
		
		Logger.printUserLog("");
	}
	
	public void execute(String[] args, String subCmdName)
	{
		Options options = getOptions();
		CommandLineParser cmdLineParser = new PosixParser();
		
		try
		{
			CommandLine cmdLine = cmdLineParser.parse(options, args, stopAtNonOption);
			CommandArguments cmdArgs = parse(cmdLine);
			cmdArgs.setOutRoot(cmdLine.getOptionValue(OPT_OUT, OPT_OUT_DEFAULT));
			
			Logger.setLogFiles(cmdArgs.getOutRoot());
			Logger.hasUserLogTag(false);
//			Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
			Logger.hasUserLogTag(true);

			printOptionsInEffect(cmdLine, subCmdName);
			
			CommandImpl cmdImpl = createCommandImpl();
			cmdImpl.preExecute();
			cmdImpl.execute(cmdArgs);
			cmdImpl.postExecute();
		}
		catch (ParseException e)
		{
			Logger.printUserError(e.getMessage());
			System.exit(1);
		}
		catch (CommandArgumentException e)
		{
			Logger.printUserError(e.getMessage());
			System.exit(1);
		}
	}
	
	protected void setStopAtNonOption(boolean stopAtNonOption)
	{
		this.stopAtNonOption = stopAtNonOption;
	}

	protected int parseIntOptionValue(CommandLine cmdLine, String opt, String defaultOptVal) throws CommandArgumentException
	{
		int value;
		try
		{
			value = Integer.parseInt(cmdLine.getOptionValue(opt, defaultOptVal));
		}
		catch (NumberFormatException e)
		{
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid integer.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}
	
	protected long parseLongOptionValue(CommandLine cmdLine, String opt) throws CommandArgumentException
	{
		return parseLongOptionValue(cmdLine, opt, null);
	}
	
	protected long parseLongOptionValue(CommandLine cmdLine, String opt, String defaultOptVal) throws CommandArgumentException
	{
		long value;
		try
		{
			value = Long.parseLong(cmdLine.getOptionValue(opt, defaultOptVal));
		}
		catch (NumberFormatException e)
		{
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid integer.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}
	
	protected long parseLongOptionValueInRange(CommandLine cmdLine, String opt, String defaultOptVal, long min, long max) throws CommandArgumentException
	{
		long value = parseLongOptionValue(cmdLine, opt, defaultOptVal);
		if (value < min || value > max)
		{
			throw new CommandArgumentException("--" + opt + " must be no smaller than " + min + " and no larger than " + max);
		}
		return value;
	}
	
	protected double parseDoubleOptionValue(CommandLine cmdLine, String opt) throws CommandArgumentException
	{
		return parseDoubleOptionValue(cmdLine, opt, null);
	}
	
	protected double parseDoubleOptionValue(CommandLine cmdLine, String opt, String defaultOptVal) throws CommandArgumentException
	{
		double value;
		try
		{
			value = Double.parseDouble(cmdLine.getOptionValue(opt, defaultOptVal));
		}
		catch (NumberFormatException e)
		{
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid floating point number.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}
	
	protected double parseDoubleOptionValueInRange(CommandLine cmdLine, String opt, String defaultOptVal, double min, double max) throws CommandArgumentException
	{
		double value = parseDoubleOptionValue(cmdLine, opt, defaultOptVal);
		if (value < min || value > max)
		{
			throw new CommandArgumentException("--" + opt + " must be no smaller than " + min + " and no larger than " + max);
		}
		return value;
	}
	
	protected String parseStringOptionValue(CommandLine cmdLine, String opt, String defaultOptVal) throws CommandArgumentException
	{
		String value = cmdLine.getOptionValue(opt, defaultOptVal);
		return value;
	}

	private boolean stopAtNonOption;
	
	protected static final String OPT_FILE_LONG = "file";
	protected static final String OPT_FILE_DESC = "Specify PLINK format .ped and .map files";
	
	protected static final String OPT_BFILE_LONG = "bfile";
	protected static final String OPT_BFILE_DESC = "Specify PLINK format .bed, .bim and .fam files";
	
	protected static final char OPT_OUT = 'o';
	protected static final String OPT_OUT_LONG = "out";
	protected static final String OPT_OUT_DEFAULT = "gear";
	protected static final String OPT_OUT_DESC = "Specify output root filename, default to '" + OPT_OUT_DEFAULT + "'";
	
	protected static final String OPT_SEED_LONG = "seed";
	protected static final String OPT_SEED_DEFAULT = "2012";
	protected static final String OPT_SEED_DESC = "Specify the seed of random number generator, default to " + OPT_SEED_DEFAULT;
}
