package gear.subcommands.simulationfam;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
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

	//dom effect
	public void setHsqDom(String hsqD) 
	{
		double hd = Double.parseDouble(hsqD);
		if(hd < 0 || hd > 0.99)
		{
			{
				Logger.printUserLog("hsq should be between 0 ~ 1.\n GEAR quitted.");
			}
		}
		if ((hd+hsq) > 1)
		{
			{
				Logger.printUserLog("The total heritability (hsq+hsq_dom) should be between 0 ~ 1.\n GEAR quitted.");
			}
		}
		hsqDom = hd;
	}

	public double getHsqDom()
	{
		return hsqDom;
	}

	public void setPlainDomEffect(double e)
	{
		polyDomEffect = e;
		isPlainDomEffect = true;
		isPolyDomEffect = false;
		isPolyDomEffectSort = false;
		isPolyDomEffectFile = false;
	}

	public boolean isPlainDomEffect()
	{
		return isPlainDomEffect;
	}

	public double getPolyDomEffect()
	{
		return polyDomEffect;
	}

	public void setPolyDomEffect()
	{
		isPlainDomEffect = false;
		isPolyDomEffect = true;
		isPolyDomEffectSort = false;
		isPolyDomEffectFile = false;
	}

	public boolean isPolyDomEffect()
	{
		return isPolyDomEffect;
	}

	public void setPolyDomEffectSort()
	{
		isPlainDomEffect = false;
		isPolyDomEffect = false;
		isPolyDomEffectSort = true;
		isPolyDomEffectFile = false;
	}

	public boolean isPolyDomEffectSort()
	{
		return isPolyDomEffectSort;
	}

	public void setPolyDomEffectFile(String f)
	{
		FileUtil.exists(f);
		polyDomEffectFile = f;

		isPlainDomEffect = false;
		isPolyDomEffect = false;
		isPolyDomEffectSort = false;
		isPolyDomEffectFile = true;
	}

	public boolean isPolyDomEffectFile()
	{
		return isPolyDomEffectFile;
	}

	public String getPolyDomEffectFile()
	{
		return polyDomEffectFile;
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

	public void setFreqFile(String ff)
	{
		FileUtil.exists(ff);
	}
	
	public boolean isFreqFile()
	{
		return isFreqFile;
	}

	public String getFreqFile()
	{
		return freqFile;
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

	public void setHsq(String h2)
	{
		hsq = Double.parseDouble(h2);
		if(hsq < 0 || hsq > 1.0)
		{
			Logger.printUserLog("hsq should be between 0~1.\n GEAR quitted.");
			System.exit(0);
		}
		isHsqFlag = true;
//		isQTLFlag = false;
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

	public void setRep(String r)
	{
		rep = Integer.parseInt(r);
		if (rep < 1)
		{
			Logger.printUserError("--rep should be greater than 0.");
			System.exit(0);
		}
	}

	public int getRep()
	{
		return rep;
	}

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
	private double hsq = 0.5;
	private boolean isHsqFlag = true;

	private double hsqDom = 0;
	private double polyDomEffect = 0.5;
	private boolean isPlainDomEffect = true;
	private boolean isPolyDomEffect = false;
	private boolean isPolyDomEffectSort = false;
	private boolean isPolyDomEffectFile = false;
	private String polyDomEffectFile = null;
	
	private String freqFile = null;
	private boolean isFreqFile = false;
	private int rep = 1;

}
