package gear.bluppca;

import gear.CommandArguments;

public class BlupPCACommandArguments extends CommandArguments
{
	public String getGRMBin()
	{
		return grmBin;
	}
	
	public void setGRMBin(String grmBin)
	{
		this.grmBin = grmBin;
	}
	
	private String grmBin;
	
	public String getGRMText()
	{
		return grmText;
	}
	
	public void setGRMText(String grmText)
	{
		this.grmText = grmText;
	}
	
	private String grmText;  // root name of the GRM files
	
	public String getGRM_GZ()
	{
		return grmGZ;
	}
	
	public void setGRM_GZ(String grmGZ)
	{
		this.grmGZ = grmGZ;
	}
	
	private String grmGZ;  // root name of the GRM files
	
	public String getGRM_ID()
	{
		return grmID;
	}
	
	public void setGRM_ID(String grmID)
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
