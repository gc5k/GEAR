package gear.subcommands.grm;

import gear.subcommands.CommandArguments;

public class GRMArguments extends CommandArguments
{
	private boolean isGZ = true;
	private boolean isVar = false;
	private boolean isDom = false;
	private double maf = 1e-5;

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
