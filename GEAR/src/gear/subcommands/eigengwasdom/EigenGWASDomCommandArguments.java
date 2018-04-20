package gear.subcommands.eigengwasdom;

import gear.subcommands.CommandArguments;

public class EigenGWASDomCommandArguments extends CommandArguments 
{
	private String pheFile;
	private int[] mPheno = {0};

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
}