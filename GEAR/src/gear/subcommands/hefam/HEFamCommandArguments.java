package gear.subcommands.hefam;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class HEFamCommandArguments extends CommandArguments
{

	public void setSD()
	{
		isSD = true;
		isSS = false;
		isCP = false;
	}
	
	public boolean isSD()
	{
		return isSD;
	}

	public void setSS()
	{
		isSD = false;
		isSS = true;
		isCP = false;
	}

	public boolean isSS()
	{
		return isSS;
	}
	
	public void setCP()
	{
		isSD = false;
		isSS = false;
		isCP = true;
	}
	
	public boolean isCP()
	{
		return isCP;
	}
	
	public String getIBDGZ()
	{
		return ibdGZ;
	}

	public void setIBDGZ(String ibdGZ)
	{
		this.ibdGZ = ibdGZ;
	}

	public String getIBDID()
	{
		return ibdID;
	}
	
	public void setIBDID(String ibdID)
	{
		this.ibdID = ibdID;
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
		this.pheIdx[0] = Integer.parseInt(pIdx);
		if (this.pheIdx[0] < 1)
		{
			Logger.printUserLog("Phenotype index should be greater than 1.\nGEAR quitted");
			System.exit(1);
		}
		this.pheIdx[0]--;
	}

	public int[] getPhenotypeIdx()
	{
		return pheIdx;
	}

	public void setIBDList(String grmL) 
	{
		ibdList = grmL;
		FileUtil.exists(ibdList);
		ArrayList<String> gFile = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(ibdList, "IBD-list");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			gFile.add(tokens[0]);
		}
		IBDList = gFile.toArray(new String[0]);
		isIBDList = true;
	}

	public boolean isIBDList()
	{
		return isIBDList;
	}
	
	public String[] getIBDList()
	{
		return IBDList;
	}
	
	public void setScale(boolean flag) 
	{
		isScale = flag;
	}

	public boolean isScale()
	{
		return isScale;
	}

	public void setCovFile(String cFile) 
	{
		FileUtil.exists(cFile);
		covFile = cFile;
	}

	public String getCovFile()
	{
		return covFile;
	}

	public void setCovNumber(String[] cIdx) 
	{
		covIdx = new int[cIdx.length];
		for (int i = 0; i < covIdx.length; i++)
		{
			covIdx[i] = Integer.parseInt(cIdx[i]);
			if (covIdx[i] < 1)
			{
				Logger.printUserLog(covIdx[i] +"< 1. Incorrect index for covar-number.");
				Logger.printUserLog("GEAR quittted.");
				System.exit(1);
			}
			covIdx[i]--;
		}
	}
	
	public int[] getCovNumber()
	{
		return covIdx;
	}

	public void setKeepFile(String kFile) 
	{
		FileUtil.exists(kFile);
		keepFile = kFile;
	}
	
	public String getKeepFile()
	{
		return keepFile;
	}

	public void setIBDcutoff(String gc) 
	{
		this.isIBDCut = true;
		this.ibdCutoff = Double.parseDouble(gc);
	}

	public double getIBDcutoff()
	{
		return ibdCutoff;
	}

	public boolean isIBDcut()
	{
		return isIBDCut;
	}

	public void setJackknife() 
	{
		isJackknife = true;
	}

	public boolean isJackknife()
	{
		return isJackknife;
	}

	private boolean isJackknife = false;
	private String covFile = null;
	private int[] covIdx = {0};
	private String ibdGZ;  // root name of the GRM files
	private String ibdID;
	private String ibdList;
	private boolean isIBDList;
	private String[] IBDList;
	private boolean isScale;

	private String pheFile = null;
	private int[] pheIdx = {0};
	private String keepFile = null;

	private boolean isSD = true;
	private boolean isSS = false;
	private boolean isCP = false;

	private boolean isIBDCut = false;
	private double ibdCutoff = 0;
}
