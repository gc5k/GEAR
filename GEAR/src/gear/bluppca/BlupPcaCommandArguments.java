package gear.bluppca;

import gear.CommandArguments;

public class BlupPcaCommandArguments extends CommandArguments
{
	public String getGrmBin()
	{
		return grmBin;
	}
	
	public void setGrmBin(String grmBin)
	{
		this.grmBin = grmBin;
	}
	
	private String grmBin;
	
	public String getGrmText()
	{
		return grmText;
	}
	
	public void setGrmText(String grmText)
	{
		this.grmText = grmText;
	}
	
	private String grmText;  // root name of the GRM files
	
	public String getGrmGZ()
	{
		return grmGZ;
	}
	
	public void setGrmGZ(String grmGZ)
	{
		this.grmGZ = grmGZ;
	}
	
	private String grmGZ;  // root name of the GRM files
	
	public String getGrmID()
	{
		return grmID;
	}
	
	public void setGrmID(String grmID)
	{
		this.grmID = grmID;
	}
	
	private String grmID;
	
	public String getPhenotypeFile()
	{
		return pheFile;
	}
	
	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}
	
	private String pheFile;
}
