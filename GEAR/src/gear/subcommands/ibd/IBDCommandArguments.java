package gear.subcommands.ibd;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class IBDCommandArguments extends CommandArguments
{
	private String chr;
	private boolean chrFlag = false;
	private boolean isGZ = true;
	private boolean scale = false;

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

	public String getChr()
	{
		return this.chr;
	}

	public boolean isChrFlagOn()
	{
		return this.chrFlag;
	}

	public void setScale() 
	{
		scale = true;
	}
	
	public boolean isScale()
	{
		return scale;
	}

}
