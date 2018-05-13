package gear.subcommands.fastpca;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class FastPCACommandArguments extends CommandArguments
{
//	public String getGrmBin()
//	{
//		return grmBin;
//	}
//
//	public void setGrmBin(String grmBin)
//	{
//		this.grmBin = grmBin;
//	}
//
//	public String getGrmText()
//	{
//		return grmText;
//	}
//	
//	public void setGrmText(String grmText)
//	{
//		this.grmText = grmText;
//	}
//	
//	private String grmText;  // root name of the GRM files
//	
//	public String getGrmGZ()
//	{
//		return grmGZ;
//	}
//	
//	public void setGrmGZ(String grmGZ)
//	{
//		this.grmGZ = grmGZ;
//	}
//	
//	private String grmGZ;  // root name of the GRM files
//	
//	public String getGrmID()
//	{
//		return grmID;
//	}
//
//	public void setGrmID(String grmID)
//	{
//		this.grmID = grmID;
//	}
	
	public void setEV(String ev)
	{
		this.ev = Integer.parseInt(ev);
	}
	
	public int getEV()
	{
		return ev;
	}

	public void setProp(String p)
	{
		prop = Double.parseDouble(p);
		if (prop < 0 & prop > 1)
		{
			Logger.printUserLog("proption should be > 0 & <=1.\n");
			Logger.printUserLog("GEAR quitted.");
			System.exit(0);
		}
	}

	public double getProp()
	{
		return prop;
	}

	public void setAdjVar()
	{
		isAdjVar = true;
	}

	public boolean isAdjVar()
	{
		return isAdjVar;
	}

	private int ev = 10;
	private double prop = 1;
	private boolean isAdjVar = false;
}
