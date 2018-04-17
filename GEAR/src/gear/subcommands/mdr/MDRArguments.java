package gear.subcommands.mdr;

import gear.subcommands.CommandArguments;

public class MDRArguments extends CommandArguments
{
	private String pheFile;
	private int mPheno;
	private int[] covIdx;
	private boolean isCC;

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

	public int getMpheno()
	{
		return this.mPheno;
	}

	public void setCovIndexes(String[] s)
	{
		covIdx = new int[s.length];
		for(int i = 1; i <covIdx.length; i++)
		{
			covIdx[i] = Integer.parseInt(s[i]);
		}
	}
	
	public int[] getCovIndexes()
	{
		return covIdx;
	}

	public void setCC(boolean hasOption)
	{
		isCC = hasOption;	
	}
	
	public int isCC()
	{
		return isCC ? 1:0; //1 for cc, 0 for qt; 
	}
}