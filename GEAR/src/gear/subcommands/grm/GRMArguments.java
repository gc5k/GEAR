package gear.subcommands.grm;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class GRMArguments extends CommandArguments
{
	private int chr;
	private boolean chrFlag = false;
	private boolean isGZ = true;
	private boolean isVar = false;

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

	public void setVar()
	{
		isVar = true;
	}

	public boolean isVar()
	{
		return isVar;
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
