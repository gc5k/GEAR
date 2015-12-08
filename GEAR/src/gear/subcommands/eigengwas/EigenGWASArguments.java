package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class EigenGWASArguments extends CommandArguments 
{
	private String pheFile;
	private int chr;
	private boolean chrFlag = false;
	private int mPheno;

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
		this.mPheno = (i - 1);
	}

	public int getMphneo()
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
}