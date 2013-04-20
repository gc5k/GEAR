package gear.util;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.logging.Level;

public class Logger
{

	private static PrintWriter userLogWriter;
	private static java.util.logging.Logger devLogger;

	static
	{
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable()
		{
			public void run()
			{
				if (userLogWriter != null)
				{
					userLogWriter.close();
					userLogWriter = null;
				}
			}
		}));

		devLogger = java.util.logging.Logger.getLogger("GEAR");
		devLogger.setUseParentHandlers(false);
		devLogger.setLevel(java.util.logging.Level.ALL);
	}

	public static void setLogFiles(String namePrefix)
	{
		String userLogFileName = namePrefix + ".log";
		try
		{
			userLogWriter = new PrintWriter(userLogFileName);
		} catch (FileNotFoundException e)
		{
			handleException(e, "Unable to create the log file '"
					+ userLogFileName + "'");
		}

		String devLogFileName = namePrefix + "_dev.log";
		try
		{
			java.util.logging.FileHandler devLogHandler = new java.util.logging.FileHandler(
					devLogFileName);
			devLogHandler.setFormatter(new java.util.logging.SimpleFormatter());
			devLogger.addHandler(devLogHandler);
		} catch (IOException e)
		{
			handleException(e, "Unable to create the log file '"
					+ devLogFileName + "'");
		}
	}

	public static void printUserLog(String msg)
	{
		System.out.println(msg);
		if (userLogWriter != null)
		{
			userLogWriter.println(msg);
		}
	}

	public static void printUserError(String msg)
	{
		System.err.println(msg);
		if (userLogWriter != null)
		{
			userLogWriter.println(msg);
		}
	}

	public static java.util.logging.Logger getDevLogger()
	{
		return devLogger;
	}

	public static void handleException(Exception e, String msg)
	{
		printUserError(msg);
		printUserError("Exception Message: " + e.getMessage());
		devLogger.log(Level.SEVERE, msg, e);
		System.exit(1);
	}

}
