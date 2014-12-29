package gear.subcommands.simulationfam;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public final class SimuFamilyCommandArguments extends CommandArguments
{
	public int getNumberOfFamilies()
	{
		return numFams;
	}
	
	public void setNumberOfFamilies(int numFams)
	{
		this.numFams = numFams;
	}
	
	public int getNumberOfMarkers()
	{
		return numMarkers;
	}
	
	public void setNumberOfMarkers(int numMarkers)
	{
		this.numMarkers = numMarkers;
	}

//start ld
	//1 plain ld
	public void setLD(double l)
	{
		this.ld = l;
		ldFlag = true;
		RandldFlag = false;
	}

	public double getLD()
	{
		return ld;
	}

	public boolean isPlainLD()
	{
		return ldFlag;
	}

	//2 unif ld
	public void setRandLD()
	{
		ldFlag = false;
		RandldFlag = true;
	}
	
	public boolean isRandLD()
	{
		return RandldFlag;
	}

//start maf
	//1 plain maf
	public void setMAF(double m)
	{
		maf = m;
		PlainmafFlag = true;
		UnifmafFlag = false;
	}

	public double getMAF()
	{
		return maf;
	}

	public boolean isPlainMAF()
	{
		return PlainmafFlag;
	}

	//2 unif maf
	public void setUnifMAF()
	{
		PlainmafFlag = false;
		UnifmafFlag = true;
	}

	public boolean isUnifMAF()
	{
		return UnifmafFlag;
	}

//start rec
	//1 plain rec
	public void setRec(double r)
	{
		this.rec = r;
		PlainrecFlag = true;
		SexrecFlag = false;
		UnifrecFlag = false;
	}
	
	public double getRec()
	{
		return rec;
	}

	public boolean isPlainRec()
	{
		return PlainrecFlag;
	}

	//2 sex rec
	public void setRecSex(String[] rs)
	{
		recSex[0] = Double.parseDouble(rs[0]);
		recSex[1] = Double.parseDouble(rs[1]);
		PlainrecFlag = false;
		SexrecFlag = true;
		UnifrecFlag = false;
	}
	
	public double[] getRecSex()
	{
		return recSex;
	}

	public boolean isSexRec()
	{
		return SexrecFlag;
	}

	//3 rand rec
	public void setRecRandFlag()
	{
		this.UnifrecFlag = true;
		this.SexrecFlag = false;
		this.PlainrecFlag = false;
	}

	public boolean isRandRec()
	{
		return UnifrecFlag;
	}

	public void setQTLFile(String f)
	{
		qtlFile = f;
		isQTLFlag = true;
		isHsqFlag = false;
	}

//start h2
	public String getQTLFile()
	{
		return qtlFile;
	}

	public boolean isQTL()
	{
		return isQTLFlag;
	}

	public void setHsq(String h2)
	{
		hsq = Double.parseDouble(h2);
		if(hsq < 0 || hsq > 1.0)
		{
			Logger.printUserLog("hsq should be between 0~1.\n GEAR quitted.");
			System.exit(0);
		}
		isHsqFlag = true;
		isQTLFlag = false;
	}

	public double getHsq()
	{
		return hsq;
	}

	public boolean isHsq()
	{
		return isHsqFlag;
	}

	
	public boolean getMakeBed()
	{
		return makeBed;
	}
	
	public void setMakeBed()
	{
		this.makeBed = true;
	}

	private int numFams = 100;
	private int numMarkers = 100;
	private boolean makeBed = false;
	
	private double maf = 0.5;
	private boolean PlainmafFlag = true;
	private boolean UnifmafFlag = false; 

	private double ld = 0;
	private boolean ldFlag = true;
	private boolean RandldFlag = false;

	private double rec = 0.5;
	private double[] recSex = {0.5, 0.5};
	private boolean PlainrecFlag = true;
	private boolean SexrecFlag = false;
	private boolean UnifrecFlag = false;

	private String qtlFile = null;
	private boolean isQTLFlag = false;
	private double hsq = 0.5;
	private boolean isHsqFlag = true;
}
