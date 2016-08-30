package gear.subcommands.oath.oathbus;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class OATHBusCommandArguments extends CommandArguments 
{
	private String covFile = null;
	private int[] covIdx = {0};

	private int chr;
	private boolean chrFlag = false;

	private String pheFile;
	private int[] mPheno;

	private String keepFile;

	private double maf = 0;
	public String getPhenotypeFile()
	{
		return this.pheFile;
	}

	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}

	public void setPhentypeIndex(String mPhe)
	{
		this.mPheno = new int[1];
		mPheno[0] = Integer.parseInt(mPhe) - 1;
		if (mPheno[0] < 0)
		{
			Logger.printUserLog("Phenotype index should be greater than 0.\nGEAR quitted.");
			System.exit(1);
		}
	}

	public int[] getMpheno()
	{
		return this.mPheno;
	}

	public void setChr(String c)
	{
		this.chr = Integer.parseInt(c);
		if (this.chr < 1)
		{
			Logger.printUserLog("Chromosome should be greater than 0.\n GEAR quitted.");
			System.exit(1);
		}
		this.chrFlag = true;
	}

	public int getChr()
	{
		return this.chr;
	}

	public boolean isChrFlagOn()
	{
		return this.chrFlag;
	}

	public void setKeeFile(String kFile) 
	{
		FileUtil.exists(kFile);
		keepFile = kFile;
	}
	
	public String getKeepFile()
	{
		return keepFile;
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
	
	public void setMAF(String mf) 
	{
		maf = Double.parseDouble(mf);
	}
	
	public double getMAF()
	{
		return maf;
	}

}
