package gear.profile;

import gear.CommandArguments;

public final class ProfileCommandArguments extends CommandArguments
{
	public String getScoreFile()
	{
		return scoreFile;
	}
	
	public void setScoreFile(String scoreFile)
	{
		this.scoreFile = scoreFile;
	}
	
	public String getQScoreFile()
	{
		return qScoreFile;
	}
	
	public void setQScoreFile(String qScoreFile)
	{
		this.qScoreFile = qScoreFile;
	}
	
	public String getQRangeFile()
	{
		return qRangeFile;
	}
	
	public void setQRangeFile(String qRangeFile)
	{
		this.qRangeFile = qRangeFile;
	}
	
	public String getFile()
	{
		return file;
	}
	
	public void setFile(String file)
	{
		this.file = file;
	}
	
	public String getBFile()
	{
		return bfile;
	}
	
	public void setBFile(String bfile)
	{
		this.bfile = bfile;
	}
	
	public String getMachDosageFile()
	{
		return machDosageFile;
	}
	
	public void setMachDosageFile(String machDosageFile)
	{
		this.machDosageFile = machDosageFile;
	}
	
	public String getMachInfoFile()
	{
		return machInfoFile;
	}
	
	public void setMachInfoFile(String machInfoFile)
	{
		this.machInfoFile = machInfoFile;
	}
	
	public String getMachDosageBatch()
	{
		return machDosageBatch;
	}
	
	public void setMachDosageBatch(String machDosageBatch)
	{
		this.machDosageBatch = machDosageBatch;
	}
	
	public String getMachInfoBatch()
	{
		return machInfoBatch;
	}
	
	public void setMachInfoBatch(String machInfoBatch)
	{
		this.machInfoBatch = machInfoBatch;
	}
	
	public CoeffModel getCoeffModel()
	{
		return coeffModel;
	}
	
	public void setCoeffModel(CoeffModel coeffModel)
	{
		this.coeffModel = coeffModel;
	}
	
	public boolean getIsSameAsPlink()
	{
		return isSameAsPlink;
	}
	
	public void setIsSameAsPlink(boolean isSameAsPlink)
	{
		this.isSameAsPlink = isSameAsPlink;
	}
	
	public String getResultFile()
	{
		return resultFile;
	}
	
	public void setResultFile(String resultFile)
	{
		this.resultFile = resultFile;
	}
	
	public boolean getIsLogit()
	{
		return isLogit;
	}
	
	public void setIsLogit(boolean isLogit)
	{
		this.isLogit = isLogit;
	}
	
	public boolean getIsAutoFlip()
	{
		return isAutoFlip;
	}
	
	public void setIsAutoFlip(boolean isAutoFlip)
	{
		this.isAutoFlip = isAutoFlip;
	}
	
	public boolean getIsWeighted()
	{
		return isWeighted;
	}
	
	public void setIsWeighted(boolean isWeighted)
	{
		this.isWeighted = isWeighted;
	}
	
	public boolean getIsKeepATGC()
	{
		return isKeepATGC;
	}
	
	public void setIsKeepATGC(boolean isKeepATGC)
	{
		this.isKeepATGC = isKeepATGC;
	}
	
	private String scoreFile;
	private String qScoreFile;
	private String qRangeFile;
	private String file;
	private String bfile;
	private String machDosageFile;
	private String machInfoFile;
	private String machDosageBatch;
	private String machInfoBatch;
	private CoeffModel coeffModel;
	private boolean isSameAsPlink;
	private String resultFile;
	private boolean isLogit;
	private boolean isAutoFlip;
	private boolean isWeighted;
	private boolean isKeepATGC;
}
