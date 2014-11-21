package gear.subcommands.weightedmeta;

import java.util.ArrayList;

import gear.gwassummary.GWASConstant;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class WeightedMetaArguments  extends CommandArguments
{
	public void setKeepCohortFile(String KFile)
	{
		FileUtil.exists(KFile);
		KeepFile = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(KFile, "Keep-cohort");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			KeepFile.add(tokens[0]);
		}

		isKeepFile = true;
		isRevFile = false;
	}

	public boolean IsKeepFile()
	{
		return isKeepFile;
	}

	public String[] getKeepFile()
	{
		return KeepFile.toArray(new String[0]);
	}
	
	public void setRemoveCohortFile(String RFile)
	{
		FileUtil.exists(RFile);
		RevFile = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(RFile, "Keep-cohort");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			RevFile.add(tokens[0]);
		}

		isKeepFile = false;
		isRevFile = true;
	}

	public boolean IsRevFile()
	{
		return isRevFile;
	}
	
	public String[] getRemoveFile()
	{
		return RevFile.toArray(new String[0]);
	}

	public void setMetaBatch(String batch)
	{
		FileUtil.exists(batch);
		md = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "MetaBatch");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			md.add(tokens[0]);
		}
	}

	public void setMetaFile(String[] m)
	{
		md = NewIt.newArrayList();
		for (int i = 0; i < m.length; i++)
		{
			FileUtil.exists(m[i]);
			md.add(m[i]);
		}
	}

	public String[] getMetaFile() 
	{
		return md.toArray(new String[0]);
	}

	public void setGZ(boolean flag)
	{
		isGZ = flag;
	}
	
	public boolean isGZ()
	{
		return isGZ;
	}

	public void setCCbatch(String ccBatch)
	{
		FileUtil.exists(ccBatch);
		ArrayList<String> s = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(ccBatch, "CC Batch");

		String[] tokens = null;
		while((tokens = reader.readTokensAtLeast(2))!=null)
		{
			s.add(tokens[0]);
			s.add(tokens[1]);
		}

		String[] cc = s.toArray(new String[0]);
		setCC(cc);
	}

	public void setCC(String[] cc)
	{
		ccSize = new double[cc.length];

		for (int i = 0; i < cc.length; i++)
		{
			ccSize[i] = Double.parseDouble(cc[i]);
			if(ccSize[i] <= 1)
			{
				Logger.printUserError("The sample size should be greater than 1.");
				System.exit(0);
			}
		}
		isQT = false;
		if ( (ccSize.length/2) != md.size())
		{
			Logger.printUserLog("The cc sample size parameters [" + ccSize.length + "] do not meet the length of the meta files [" + md.size()+"].");
			System.exit(0);
		}

	}
	
	public double[] getCCsize()
	{
		return ccSize;
	}

	public void setQTbatch(String qtBatch)
	{
		ArrayList<String> s = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(qtBatch, "QT Batch");

		String[] tokens = null;
		while((tokens = reader.readTokensAtLeast(1))!=null)
		{
			s.add(tokens[0]);
		}

		String[] qt = s.toArray(new String[0]);
		setQT(qt);
	}

	public void setQT(String[] qt)
	{
		qtSize = new double[qt.length];
		
		for (int i = 0; i < qtSize.length; i++)
		{
			qtSize[i] = Double.parseDouble(qt[i]);
			if (qtSize[i] <= 1)
			{
				Logger.printUserError("The sample size should be greater than 1.");
				System.exit(0);
			}
		}
		isQT = true;
		
		if ( qtSize.length != md.size())
		{
			Logger.printUserLog("The qt sample size parameters [" + qtSize.length + "] do not meet the length of the meta files [" + md.size()+"].");
			System.exit(0);
		}
	}

	public double[] getQTsize()
	{
		return qtSize;
	}

	public boolean isQT() 
	{
		return isQT;
	}

	public void setGC()
	{
		isGC = true;
	}

	public boolean getGC()
	{
		return isGC;
	}

	public void setGCInflationOnly()
	{
		isGCInflationOnly = true;
	}

	public boolean getGCInflationOnly()
	{
		return isGCInflationOnly;
	}

	public void setATGC()
	{
		isKeepATGC = true;
	}

	public boolean isKeepATGC()
	{
		return isKeepATGC;
	}

	public String[] getKeys()
	{
		return field;
	}

	public void setKey(String[] k)
	{
		field[GWASConstant.SNP] = k[0];
		if (isQT)
		{
			field[GWASConstant.BETA] = k[1];			
		}
		else
		{
			field[GWASConstant.OR] = k[1];			
		}
		field[GWASConstant.SE] = k[2];
		field[GWASConstant.A1] = k[3];
		field[GWASConstant.A2] = k[4];

		if(k.length >5)
		{
			field[GWASConstant.CHR] = k[5];
		}
		if(k.length >6)
		{
			field[GWASConstant.BP] = k[6];
		}
		if(k.length >7)
		{
			field[GWASConstant.P] = k[7];
		}
	}

	public void setVerbose()
	{
		isVerbose = true;
	}

	public void setVerboseGZ()
	{
		isVerbose = true;
		isVerboseGZ = true;
	}

	public boolean isVerbose()
	{
		return isVerbose;
	}

	public boolean isVerboseGZ()
	{
		return isVerboseGZ;
	}

	public void setCM(String cF)
	{
		FileUtil.exists(cF);
		cmFile = cF;
		isCM = true;
	}

	public String getCMFile()
	{
		return cmFile;
	}

	public boolean isCM()
	{
		return isCM;
	}

	public void setFullSNPOnly()
	{
		isFullSNPOnly = true;
	}
	
	public boolean isFullSNPOnly()
	{
		return isFullSNPOnly;
	}

	private ArrayList<String> md;
	private boolean isGZ = false;
	private boolean isQT = true;
	private boolean isGC = false;
	private boolean isGCInflationOnly = false;
	private boolean isKeepATGC = false;
	private boolean isVerbose = false;
	private boolean isVerboseGZ = false;

	private double[] qtSize;
	private double[] ccSize;
	private String[] field = {"snp", "chr", "bp", "beta", "or", "se", "p", "a1", "a2"};
	
	private String cmFile = null;
	private boolean isCM = false;
	private boolean isFullSNPOnly = false;
	
	private ArrayList<String> KeepFile = null;
	private boolean isKeepFile = false;
	private ArrayList<String> RevFile = null;
	private boolean isRevFile = false;
}
