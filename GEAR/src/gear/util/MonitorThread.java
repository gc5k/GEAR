package gear.util;

import gear.ConstValues;

public class MonitorThread extends Thread
{
	public MonitorThread()
	{
		setPriority(MIN_PRIORITY);
	}
	
	public void run() {
		isRun = true;
		do
		{
			long totalMem = Runtime.getRuntime().totalMemory();
			peakMem = Math.max(peakMem, totalMem);
			try {
				Thread.sleep(500);
			}
			catch (InterruptedException e) {
				// Do nothing
			}
		} while (isRun);
	}
	
	public void stopMonitoring()
	{
		isRun = false;
	}
	
	public long getPeakMemory()
	{
		return peakMem;
	}
	
	public String getPeakMemoryFormatString()
	{
		if (getPeakMemory() >= ConstValues.GIGABYTE)
		{
			return String.format("%.1f GB", ((double)getPeakMemory()) / ConstValues.GIGABYTE);
		}
		else if (getPeakMemory() >= ConstValues.MEGABYTE)
		{
			return String.format("%.1f MB", ((double)getPeakMemory()) / ConstValues.MEGABYTE);
		}
		else if (getPeakMemory() >= ConstValues.KILOBYTE)
		{
			return String.format("%.1f KB", ((double)getPeakMemory()) / ConstValues.KILOBYTE);
		}
		return String.format("%d byte(s)", getPeakMemory());
	}
	
	private boolean isRun;
	private long peakMem;
}
