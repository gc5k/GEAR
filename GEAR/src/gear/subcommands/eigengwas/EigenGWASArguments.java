package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class EigenGWASArguments extends CommandArguments 
{
	private String pheFile;
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
		this.mPheno[0] = (i - 1);
	}

	public int[] getMpheno()
	{
		return this.mPheno;
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
