package gear.subcommands.profile;

import gear.subcommands.CommandArguments;

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
	
	public boolean getHasScoreHeader()
	{
		return hasScoreHeader;
	}

	public void setHasScoreHeader(boolean hasScoreHeader)
	{
		this.hasScoreHeader = hasScoreHeader;
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
	
	public CoeffModelType getCoeffModelType()
	{
		return coeffModelType;
	}
	
	public void setCoeffModelType(CoeffModelType coeffModelType)
	{
		this.coeffModelType = coeffModelType;
	}
	
	public String getCoeffModelFile()
	{
		return coeffModelFile;
	}
	
	public void setCoeffModelFile(String coeffModelFile)
	{
		this.coeffModelFile = coeffModelFile;
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
	private boolean hasScoreHeader;
	private String qScoreFile;
	private String qRangeFile;
	private String machDosageFile;
	private String machInfoFile;
	private String machDosageBatch;
	private String machInfoBatch;
	private CoeffModelType coeffModelType;
	private String coeffModelFile;
	private boolean isSameAsPlink;
	private String resultFile;
	private boolean isLogit;
	private boolean isAutoFlip;
	private boolean isWeighted;
	private boolean isKeepATGC;
}
