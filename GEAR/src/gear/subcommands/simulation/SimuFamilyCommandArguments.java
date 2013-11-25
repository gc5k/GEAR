package gear.subcommands.simulation;

import gear.subcommands.CommandArguments;

public final class SimuFamilyCommandArguments extends CommandArguments
{
	public int getNumberOfFamilies()
	{
		return numFams;
	}
	
	public void setNumberOfFamilies(int numFams)
	{
		this.numFams = numFams;
	}
	
	public int getNumberOfMarkers()
	{
		return numMarkers;
	}
	
	public void setNumberOfMarkers(int numMarkers)
	{
		this.numMarkers = numMarkers;
	}
	
	public Long getSeed()
	{
		return seed;
	}
	
	public void setSeed(Long seed)
	{
		this.seed = seed;
	}
	
	public boolean getMakeBed()
	{
		return makeBed;
	}
	
	public void setMakeBed(boolean makeBed)
	{
		this.makeBed = makeBed;
	}
	
	private int numFams;
	private int numMarkers;
	private Long seed = null;
	private boolean makeBed;
}
