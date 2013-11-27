package gear.subcommands;

import gear.util.Logger;
import gear.util.MonitorThread;

import java.util.Calendar;

public abstract class CommandImpl
{
	public void preExecute()
	{
		Logger.printUserLog("Analysis started: " + Calendar.getInstance().getTime());
		Logger.printUserLog("");
		monitor = new MonitorThread();
		monitor.start();
	}
	
	public abstract void execute(CommandArguments cmdArgs);
	
	public void postExecute()
	{
		monitor.stopMonitoring();
		Logger.printUserLog("");
		Logger.printUserLog("Analysis finished: " + Calendar.getInstance().getTime());
		Logger.printUserLog("Peak memory consumption: " + monitor.getPeakMemoryFormatString());
	}
	
	private MonitorThread monitor;
}
