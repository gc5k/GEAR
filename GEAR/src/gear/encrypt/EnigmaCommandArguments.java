package gear.encrypt;

import gear.CommandArguments;

public class EnigmaCommandArguments extends CommandArguments
{
	public long getSeed()
	{
		return seed;
	}
	
	public void setSeed(long seed)
	{
		this.seed = seed;
	}
	
	private long seed;
	
	public int getNumberOfColumns()
	{
		return numCols;
	}
	
	public void setNumberOfColumns(int numCols)
	{
		this.numCols = numCols;
	}
	
	private int numCols;
	
	public String getMapFile()
	{
		return mapFile;
	}
	
	public void setMapFile(String mapFile)
	{
		this.mapFile = mapFile;
	}
	
	private String mapFile;
}
