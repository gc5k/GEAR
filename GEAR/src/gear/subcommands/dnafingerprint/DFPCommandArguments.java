package gear.subcommands.dnafingerprint;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class DFPCommandArguments extends CommandArguments
{

	public void setBFile2(String bFile2) 
	{
		this.bFile2 = bFile2;
		FileUtil.exists(new String(this.bFile2 + ".bed"));
		FileUtil.exists(new String(this.bFile2 + ".bim"));
		FileUtil.exists(new String(this.bFile2 + ".fam"));
	}
	
	public String getBFile2()
	{
		return bFile2;
	}

	public void setLowCutoff(double lowCutoff)
	{
		this.lowCutoff = lowCutoff;
	}
	
	public double getLowCutoff()
	{
		return lowCutoff;
	}
	
	public void setHighCutoff(double highCutoff)
	{
		this.highCutoff = highCutoff;
	}
	
	public double getHighCutoff()
	{
		return highCutoff;
	}

	public void setNumMarker(long numMarker)
	{
		this.numMarker = numMarker;
	}
	
	public long getNumMarker()
	{
		return numMarker;
	}

	private String bFile2 = null;

	private double lowCutoff;
	private double highCutoff;
	private long numMarker;
}
