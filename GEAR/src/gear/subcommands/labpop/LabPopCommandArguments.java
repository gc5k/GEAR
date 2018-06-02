package gear.subcommands.labpop;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math.random.RandomDataImpl;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class LabPopCommandArguments extends CommandArguments
{
	public int getSampleSize()
	{
		return sampleSize;
	}

	public void setSampleSize(int ss)
	{
		this.sampleSize = ss;
	}

	public int getNumberOfMarkers()
	{
		return numMarkers;
	}

	public void setNumberOfMarkers(int numMarkers)
	{
		this.numMarkers = numMarkers;
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
	
	public int getNullMarkerNum()
	{
		return nullM;
	}

// start rec
	// 1 plain rec
	public void setRec(double r)
	{
		rec = new double[getNumberOfMarkers()];
		Arrays.fill(rec, r);
		rec[0] = 0.5;
	}

	public double[] getRec()
	{
		return rec;
	}

	// 2 rand rec
	public void setRecRand() 
	{
		rec = new double[getNumberOfMarkers()];
		rnd.reSeed(getSeed());
		for(int i = 0; i < rec.length; i++)
		{
			rec[i] = rnd.nextUniform(0.01, 0.5);
		}
		rec[0] = 0.5;
	}

	// 3 rec file
	public void setRecFile(String f)
	{
		FileUtil.exists(f);
		RecFile = f;
		
		BufferedReader reader = FileUtil.FileOpen(f);
		rec = new double[getNumberOfMarkers()];

		int c = 0;
		String line = null;
		try
		{
			while ((line = reader.readLine()) != null)
			{
				if(c >= getNumberOfMarkers())
				{
					Logger.printUserLog("Have already read " + getNumberOfMarkers() + " recombination fractions.  Ignore the rest of the content in '" + f + "'.");
					break;
				}

				line.trim();
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 1) continue;
				rec[c++] = Double.parseDouble(l[0]);
				if (rec[c] > 0.5 || rec[c] < 0)
				{
					Logger.printUserLog("incorrect recombination fraction : '" + rec[c] + "' in line " + c + ".\n Gear quit.");
					System.exit(0);
				}
			}
			reader.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e,
						"An exception occurred when reading the recombination fraction file '"
								+ f + "'.");
		}

		if (rec.length < getNumberOfMarkers())
		{
			Logger.printUserLog("The number of recombination fraction is fewer than the number of markers.\nGEAR quit.");
			System.exit(0);
		}
		rec[0] = 0.5;
	}
	
	public String getRecFile()
	{
		return RecFile;
	}

//heritability	
	public void setHsq(double h2)
	{
		hsq = h2;
//		isHsqFlag = true;
	}

	public double getHsq()
	{
		return hsq;
	}

	public void setHsqDom(String hsqD) 
	{
		double hd = Double.parseDouble(hsqD);
		if (hd < 0 || hd > 0.99)
		{
			{
				Logger.printUserLog("hsq should be between 0 ~ 1.\n GEAR quit.");
			}
		}
		if ((hd+hsq) > 1)
		{
			{
				Logger.printUserLog("The total heritability (hsq+hsq_dom) should be between 0 ~ 1.\n GEAR quit.");
			}
		}
		hsqDom = hd;
	}

	public double getHsqDom()
	{
		return hsqDom;
	}

//	public boolean isHsq()
//	{
//		return isHsqFlag;
//	}

//polygenic effects

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
	
//dom effect

	public void setPolyDomEffectFile(String f)
	{
		FileUtil.exists(f);
		polyDomEffectFile = f;
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


//bed	
	public boolean getMakeBed()
	{
		return makeBed;
	}
	
	public void setMakeBed()
	{
		this.makeBed = true;
	}

//population
	public void setBC()
	{
		isBC = true;
		isF2 = false;
		isDH = false;
		isRIL = false;
		isIF2 = false;
	}

	public boolean isBC()
	{
		return isBC;
	}
	
	public void setF2()
	{
		isBC = false;
		isF2 = true;
		isDH = false;
		isRIL = false;
		isIF2 = false;
	}
	public boolean isF2()
	{
		return isF2;
	}
	
	public void setDH()
	{
		isBC = false;
		isF2 = false;
		isDH = true;
		isRIL = false;
		isIF2 = false;
	}

	public boolean isDH()
	{
		return isDH;
	}

	public void setRIL()
	{
		isBC = false;
		isF2 = false;
		isDH = false;
		isRIL = true;
		isIF2 = false;
	}
	public boolean isRIL()
	{
		return isRIL;
	}

	public void setIF2()
	{
		isBC = false;
		isF2 = false;
		isDH = false;
		isRIL = false;
		isIF2 = true;
	}

	public boolean isIF2()
	{
		return isIF2;
	}

	public void setReplication(String rep)
	{
		replication = Integer.parseInt(rep);
		if (replication < 1)
		{
			Logger.printUserLog("Replication should be greater than 0. Gear quit.");
			System.exit(0);
		}
	}
	
	public int getReplication()
	{
		return replication;
	}

	public void set1234mode()
	{
		is1234mode = true;
		isATGCmode = false;
	}
	public boolean is1234mode()
	{
		return is1234mode;
	}
	
	public void setATGCmode()
	{
		isATGCmode = true;
		is1234mode = false;
	}
	public boolean isATGCmode()
	{
		return isATGCmode;
	}

	public void setHSQB(double e)
	{
		hsqB = e;
	}

	public double getHSQB()
	{
		return hsqB;
	}

	public void setVerbose()
	{
		isVerbose = true;
	}
	
	public boolean isVerbose()
	{
		return isVerbose;
	}

	private boolean isVerbose = false;
	private boolean isBC = false;
	private boolean isF2 = true;
	private boolean isDH = false;
	private boolean isRIL = false;
	private boolean isIF2 = false;

	private boolean isPolyEffectFile = false;
	private String polyEffectFile = null;

	private boolean isPolyDomEffectFile = false;
	private String polyDomEffectFile = null;

	private int sampleSize = 100;
	private int numMarkers = 100;
	private int nullM = 0;
	private boolean makeBed = false;

	private double[] rec;
	private String RecFile = null;

	private double hsq = 0.5;
	private double hsqDom = 0;
	private double hsqB = -1;

//	private boolean isHsqFlag = true;
	
	private int replication = 1;

	private RandomDataImpl rnd = new RandomDataImpl();		

	private boolean is1234mode = false;
	private boolean isATGCmode = true;
}
