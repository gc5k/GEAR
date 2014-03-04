package gear.subcommands.metawatchdog.decrypt;

import gear.subcommands.CommandArguments;

public class MetaWatchdog2CommandArguments extends CommandArguments
{
	public String getDataset1()
	{
		return dataset1;
	}
	public void setDataset1(String dataset1)
	{
		this.dataset1 = dataset1;
	}
	public String getDataset2()
	{
		return dataset2;
	}
	public void setDataset2(String dataset2)
	{
		this.dataset2 = dataset2;
	}
	public void setSquare(float sd)
	{
		squareFlag = true;
		squareDis = sd;
	}
	public boolean getSquare()
	{
		return squareFlag;
	}
	public float getSquareDis()
	{
		return squareDis;
	}
	public float getCutoff()
	{
		return cutoff;
	}
	public void setCutoff(float cutoff)
	{
		this.cutoff = cutoff;
	}
	private String dataset1;
	private String dataset2;
	private float cutoff;
	private boolean squareFlag = false;
	private float squareDis;
}
