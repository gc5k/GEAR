package gear;

public abstract class CommandArguments
{
	public void setFile(String file)
	{
		this.file = file;
	}
	
	public String getFile()
	{
		return file;
	}
	
	public String getPed()
	{
		return getFile() + ".ped";
	}
	
	public String getMap()
	{
		return getFile() + ".map";
	}
	
	private String file;
	
	public String getBFile()
	{
		return bfile;
	}
	
	public void setBFile(String bfile)
	{
		this.bfile = bfile;
	}
	
	public String getBed()
	{
		return getBFile() + ".bed";
	}
	
	public String getBim()
	{
		return getBFile() + ".bim";
	}
	
	public String getFam()
	{
		return getBFile() + ".fam";
	}
	
	private String bfile;
	
	public String getOutRoot()
	{
		return outRoot;
	}
	
	public void setOutRoot(String outRoot)
	{
		this.outRoot = outRoot;
	}
	
	private String outRoot;
}
