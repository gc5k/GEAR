package gear.help;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import gear.Command;
import gear.CommandArguments;
import gear.CommandImpl;
import gear.Gear;
import gear.util.Logger;

public class HelpCommandImpl extends CommandImpl
{
	@Override
	public void preExecute()
	{
		Logger.hasUserLogTag(false);
	}
	
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		HelpCommandArguments helpCmdArgs = (HelpCommandArguments)cmdArgs;
		String[] subcmds = helpCmdArgs.getSubcommands();
		
		if (subcmds.length == 0)
		{
			printTopLevelHelp();
		}
		else
		{
			for (String subcmd : subcmds)
			{
				printSubcommandHelp(subcmd);
			}
		}
	}
	
	@Override
	public void postExecute()
	{
	}

	private void printTopLevelHelp()
	{
		Logger.printUserLog("Usage: gear <subcommand> [options] [args]");
		Logger.printUserLog("Type 'gear help <subcommand>' for help on a specific subcommand.");
		Logger.printUserLog("");
		Logger.printUserLog("Available subcommands:");
		printAvailableSubcommands();
		Logger.printUserLog("");
	}
	
	private void printAvailableSubcommands()
	{
		Set<Command> cmdSet = Gear.INSTANCE.getCommandSet();
		Iterator<Command> cmdIter = cmdSet.iterator();
		String[] names = new String[cmdSet.size()];
		String[] descs = new String[cmdSet.size()];
		int longestNameEndPos = 0;
		
		for (int cmdIdx = 0; cmdIter.hasNext(); ++cmdIdx)
		{
			Command cmd = cmdIter.next();
			names[cmdIdx] = "   " + getNameAndAliasesString(cmd);
			
			if (names[cmdIdx].length() > longestNameEndPos)
			{
				longestNameEndPos = names[cmdIdx].length();
			}
			
			descs[cmdIdx] = cmd.getDescription();
		}
		
		int descStartPos = longestNameEndPos + 3;
		
		for (int cmdIdx = 0; cmdIdx < cmdSet.size(); ++cmdIdx)
		{
			String cmdEntry = names[cmdIdx];
			for (int spaceCnt = 0; spaceCnt < descStartPos - names[cmdIdx].length(); ++spaceCnt)
			{
				cmdEntry += " ";
			}
			cmdEntry += descs[cmdIdx];
			Logger.printUserLog(cmdEntry);
		}
	}
	
	private void printSubcommandHelp(String sNameOrAlias)
	{
		Command cmd = Gear.INSTANCE.getCommand(sNameOrAlias);
		if (cmd == null)
		{
			Logger.printUserLog("\"" + sNameOrAlias + "\": unknown subcommand");
			Logger.printUserLog("");
		}
		else
		{
			Logger.printUserLog(getNameAndAliasesString(cmd) + ": " + cmd.getLongDescription());
			Logger.printUserLog("");
			
			if (cmd.getFullDescription().length() != 0)
			{
				Logger.printUserLog(cmd.getFullDescription());
				Logger.printUserLog("");
			}
			
			Options options = cmd.getOptions();
			if (!options.getOptions().isEmpty())
			{
				Logger.printUserLog("Available Options:");
				HelpFormatter helpFormatter = new HelpFormatter();
				StringWriter stringWriter = new StringWriter();
				PrintWriter printWriter = new PrintWriter(stringWriter);
				helpFormatter.printOptions(printWriter,
				                           helpFormatter.getWidth(),
				                           options,
				                           helpFormatter.getLeftPadding(),
				                           helpFormatter.getDescPadding());
				Logger.printUserLog(stringWriter.toString());
				Logger.printUserLog("");
			}
		}
	}
	
	private String getNameAndAliasesString(Command cmd)
	{
		String result = cmd.getName();
		
		Set<String> aliases = cmd.getAliases();
		if (!aliases.isEmpty())
		{
			result += " (";
			Iterator<String> aliasIter = aliases.iterator();
			boolean firstAlias = true;
			while (aliasIter.hasNext())
			{
				if (firstAlias)
				{
					firstAlias = false;
				}
				else
				{
					result += ", ";
				}
				result += aliasIter.next();
			}
			result += ")";
		}
		
		return result;
	}
}
