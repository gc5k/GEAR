package gear.subcommands.eigengwasdom;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class EigenGWASDomCommandArguments extends CommandArguments 
{
	private String pheFile;
	private int chr;
	private boolean chrFlag = false;
	private int[] mPheno = {0};
	private String keepFile = null;

	public String getPhenotypeFile()
	{
		return this.pheFile;
	}

	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}

	public void setPhentypeIndex(int i)
	{
		this.mPheno[0] = i - 1;
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
			Logger.printUserLog("Chromosome should be greater than 0.\n GEAR quitted");
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
	
	public void setKeepFile(String kFile) 
	{
		FileUtil.exists(kFile);
		keepFile = kFile;
	}
	
	public String getKeepFile()
	{
		return keepFile;
	}

}