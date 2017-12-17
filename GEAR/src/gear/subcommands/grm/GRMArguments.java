package gear.subcommands.grm;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class GRMArguments extends CommandArguments
{
	private String chr;
	private boolean chrFlag = false;
	private boolean isGZ = true;
	private boolean isVar = false;
	private boolean isDom = false;
	private double maf = 1e-5;

	public void setChr(String c)
	{
		this.chr = c;
		int chr1 = Integer.parseInt(c);
		if (chr1 < 1)
		{
			Logger.printUserLog("Chromosome should be greater than 0.\n GEAR quitted");
			System.exit(1);
		}
		this.chrFlag = true;
	}

	public void setGZ()
	{
		isGZ = true;
	}

	public void setTxt()
	{
		isGZ = false;
	}
	
	public boolean isGZ()
	{
		return isGZ;
	}

	public void setAdjVar()
	{
		isVar = true;
	}

	public boolean isAdjVar()
	{
		return isVar;
	}

	public String getChr()
	{
		return this.chr;
	}

	public boolean isChrFlagOn()
	{
		return this.chrFlag;
	}

	public void setDom() 
	{
		isDom = true;
	}

	public boolean isDom()
	{
		return isDom;
	}

	public void setMAF(String m) 
	{
		maf = Double.parseDouble(m);
	}
	
	public double getMAF()
	{
		return maf;
	}
}
