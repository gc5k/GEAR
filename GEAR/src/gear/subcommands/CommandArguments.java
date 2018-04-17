package gear.subcommands;

import java.util.HashSet;

import gear.util.NewIt;

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
	
	public String getOutRoot()
	{
		return outRoot;
	}
	
	public void setOutRoot(String outRoot)
	{
		this.outRoot = outRoot;
	}
	
	public long getSeed()
	{
		return seed;
	}

	public void setSeed(long seed)
	{
		this.seed = seed;
	}

	public void setKeepFile(String kF)
	{
		this.keepFile = kF;
	}

	public String getKeepFile()
	{
		return keepFile;
	}

	public void setRemoveFile(String rF)
	{
		this.removeFile = rF;
	}

	public String getRemoveFile()
	{
		return removeFile;
	}

	public void setExtractFile(String eF)
	{
		this.extractFile = eF;
	}

	public String getExtractFile()
	{
		return extractFile;
	}

	public void setExcludeFile(String exF)
	{
		this.excludeFile = exF;
	}

	public String getExcludeFile()
	{
		return excludeFile;
	}

	public void setChr(String[] c)
	{
		chr = NewIt.newHashSet();
		for(int i = 0; i < c.length; i++)
		{
			chr.add(c[i]);
		}
	}
	
	public HashSet<String> getChr()
	{
		return chr;
	}

	public void setNotChr(String[] c)
	{
		not_chr = NewIt.newHashSet();
		for(int i = 0; i < c.length; i++)
		{
			not_chr.add(c[i]);
		}
	}

	public HashSet<String> getNotChr()
	{
		return not_chr;
	}

	private String file;
	private String bfile;
	private String outRoot;
	private String keepFile = null;
	private String removeFile = null;
	private String extractFile = null;
	private String excludeFile = null;
	private HashSet<String> chr = null;
	private HashSet<String> not_chr = null;
	private long seed;
}
