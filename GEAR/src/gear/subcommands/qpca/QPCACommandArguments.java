package gear.subcommands.qpca;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class QPCACommandArguments extends CommandArguments
{
	public String getGrmBin()
	{
		return grmBin;
	}

	public void setGrmBin(String grmBin)
	{
		FileUtil.exists(grmBin);
		this.grmBin = grmBin;
	}

	public String getGrmText()
	{
		return grmText;
	}

	public void setGrmText(String grmText)
	{
		FileUtil.exists(grmText);
		this.grmText = grmText;
	}

	private String grmText;  // root name of the GRM files
	
	public String getGrmGZ()
	{
		return grmGZ;
	}
	
	public void setGrmGZ(String grmGZ)
	{
		FileUtil.exists(grmGZ);
		this.grmGZ = grmGZ;
	}
	
	private String grmGZ;  // root name of the GRM files
	
	public String getGrmID()
	{
		return grmID;
	}

	public void setGrmID(String grmID)
	{
		FileUtil.exists(grmID);
		this.grmID = grmID;
	}
	
	public void setEV(String ev)
	{
		this.ev = Integer.parseInt(ev);
	}
	
	public int getEV()
	{
		return ev;
	}

	private String grmBin;
	private String grmID;
	private int ev = 10;
}
