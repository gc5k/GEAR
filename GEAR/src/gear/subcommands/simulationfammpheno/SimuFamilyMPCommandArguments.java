package gear.subcommands.simulationfammpheno;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public final class SimuFamilyMPCommandArguments extends CommandArguments
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

	public int getNullMarker()
	{
		return nullM;
	}

	public void setNullMarkerNum(String nm)
	{
		nullM = Integer.parseInt(nm);
		if (nullM < 0)
		{
			Logger.printUserLog("Null marker number " + nullM + " is negative. It is set to zero.");
			nullM = 0;
		}
		if (nullM >= numMarkers)
		{
			Logger.printUserLog("Null marker number " + nullM + " should less than the number of merkers (" + numMarkers +").\n GEAR quittte.");
			System.exit(0);
		}
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

//	public void setQTLFile(String f)
//	{
//		qtlFile = f;
//		isQTLFlag = true;
//		isHsqFlag = false;
//	}

//start h2
//	public String getQTLFile()
//	{
//		return qtlFile;
//	}
//
//	public boolean isQTL()
//	{
//		return isQTLFlag;
//	}

	public void setHsq(String[] h2)
	{
		hsq = new double[h2.length];
		for(int i = 0; i < h2.length; i++)
		{
			hsq[i] = Double.parseDouble(h2[i]);
			if(hsq[i] < 0 || hsq[i] > 1.0)
			{
				Logger.printUserLog("hsq should be between 0~1. GEAR quit.");
				System.exit(0);
			}
			isHsqFlag = true;
//			isQTLFlag = false;			
		}
	}

	public double[] getHsq()
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
	
	public void setPolyEffectFile(String f)
	{
		FileUtil.exists(f);
		polyEffectFile = f;
		isPolyEffectFile = true;
	}

	public boolean isPolyEffectFile()
	{
		return isPolyEffectFile;
	}

	public String getPolyEffectFile()
	{
		return polyEffectFile;
	}

	public void setCM(String cF)
	{
		FileUtil.exists(cF);
		cmFile = cF;
	}

	public String getCMFile()
	{
		return cmFile;
	}

	public void setCME(String ceF)
	{
		FileUtil.exists(ceF);
		cmeFile = ceF;
	}

	public String getCMEFile()
	{
		return cmeFile;
	}

	private String cmFile = null;
	private String cmeFile = null;

	private boolean isPolyEffectFile = false;
	private String polyEffectFile = null;

	private int numFams = 100;
	private int numMarkers = 100;
	private int nullM = 0;

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

//	private String qtlFile = null;
//	private boolean isQTLFlag = false;
	private double[] hsq = {0.5};
	private boolean isHsqFlag = true;
}
