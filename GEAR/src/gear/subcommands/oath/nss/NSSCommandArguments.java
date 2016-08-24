package gear.subcommands.oath.nss;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class NSSCommandArguments extends CommandArguments 
{
	private String pheFile;
	private int chr;
	private boolean chrFlag = false;
	private int[] mPheno;
	private String keepFile;

	public String getPhenotypeFile()
	{
		return this.pheFile;
	}

	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}

	public void setPhentypeIndex(String[] mPhe)
	{
		this.mPheno = new int[mPhe.length];
		for (int i = 0; i < mPhe.length; i++)
		{
			mPheno[i] = Integer.parseInt(mPhe[i]) - 1;
			if (mPheno[i] < 0)
			{
				Logger.printUserLog("Phenotype index should be greater than 0.\nGEAR quitted.");
				System.exit(1);
			}
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
}
