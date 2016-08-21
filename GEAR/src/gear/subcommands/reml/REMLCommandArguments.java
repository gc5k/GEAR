package gear.subcommands.reml;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class REMLCommandArguments extends CommandArguments 
{
	public String getGrmBin()
	{
		return grmBin;
	}
	
	public void setGrmBin(String grmBin)
	{
		this.grmBin = grmBin;
	}

	public String getGrmText()
	{
		return grmText;
	}
	
	public void setGrmText(String grmText)
	{
		this.grmText = grmText;
	}
	
	public String getGrmGZ()
	{
		return grmGZ;
	}
	
	public void setGrmGZ(String grmGZ)
	{
		this.grmGZ = grmGZ;
	}
	
	public String getGrmID()
	{
		return grmID;
	}
	
	public void setGrmID(String grmID)
	{
		this.grmID = grmID;
	}
	

	public String getPhenotypeFile()
	{
		return pheFile;
	}

	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}

	public void setPhenotypeIdx(String pIdx)
	{
		this.pheIdx = Integer.parseInt(pIdx);
		if (pheIdx < 1)
		{
			Logger.printUserLog("Phenotype index should be greater than 1.\nGEAR quitted");
			System.exit(1);
		}
		this.pheIdx--;
	}

	public int getPhenotypeIdx()
	{
		return pheIdx;
	}

	public void setGrmList(String grmL) 
	{
		grmList = grmL;
		FileUtil.exists(grmList);
		ArrayList<String> gFile = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(grmList, "GRM-list");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			gFile.add(tokens[0]);
		}
		grmGZList = gFile.toArray(new String[0]);
		isGRMList = true;
	}

	public boolean isGRMList()
	{
		return isGRMList;
	}
	
	public String[] getGrmList()
	{
		return grmGZList;
	}
	
	public void setMINQUE(boolean flag) 
	{
		isMINQUE = flag;
	}
	
	public boolean isMINQUE()
	{
		return isMINQUE;
	}

	private String grmBin;
	private String grmText;  // root name of the GRM files
	private String grmGZ;  // root name of the GRM files
	private String grmID;
	private String grmList;
	private boolean isGRMList;
	private String[] grmGZList;
	private boolean isMINQUE;

	private String pheFile;
	private int pheIdx = 0;
}
