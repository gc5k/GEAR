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

	public boolean getRecRand()
	{
		return recRandFlag;
	}

	public void setRecRand(boolean f)
	{
		this.recRandFlag = f;
	}


	public double getLD()
	{
		return ld;
	}

	public void setLD(double l)
	{
		this.ld = l;
	}

	public double getRec()
	{
		return rec;
	}

	public void setRec(double r)
	{
		this.rec = r;
	}
	
	private int numFams;
	private int numMarkers;
	private Long seed = null;
	private boolean makeBed;
	private double ld = 0;
	private double rec = 0.5;
	private boolean recRandFlag = false;
}
