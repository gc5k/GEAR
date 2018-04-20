package gear.subcommands.eigengwasepi;

import gear.subcommands.CommandArguments;

public class EigenGWASEpiCommandArguments extends CommandArguments 
{
	private String pheFile;
	private int[] mPheno = {0};
	private boolean inbred = false;

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

	public void setInbred()
	{
		inbred = true;
	}

	public boolean isInbred()
	{
		return inbred;
	}
}
