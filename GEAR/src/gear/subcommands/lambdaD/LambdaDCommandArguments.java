package gear.subcommands.lambdaD;

import java.util.ArrayList;

import gear.gwassummary.GWASConstant;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class LambdaDCommandArguments extends CommandArguments
{

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

	public String[] getKeys()
	{
		return field;
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

	public double getMe()
	{
		return Me;
	}

	public void setMe(String me)
	{
		Me = Double.parseDouble(me);
	}

	public void setChr(String chr)
	{
		this.chr = Integer.parseInt(chr); 
		chrFlag = true;
	}

	public int getChr()
	{
		return chr;
	}
	
	public boolean isChr()
	{
		return chrFlag;
	}

	public boolean isBeta()
	{
		return isBeta;
	}

	public void setFrq()
	{
		isBeta = false;
		isFrq = true;
		mode = FRQ;
	}

	public boolean isFrq()
	{
		return isFrq;
	}

	public void setFst()
	{
		isBeta = false;
		isFrq = true;
		mode = FRQ;
	}

	public void setClean()
	{
		isClean = true;
	}
	
	public boolean isClean()
	{
		return isClean;
	}

	public void setRapid()
	{
		isRapid = true;
	}
	
	public boolean isRapid()
	{
		return isRapid;
	}

	public void setRandom()
	{
		isRandom = true;
	}
	
	public boolean isRandom()
	{
		return isRandom;
	}

	public void setTrim(double tr)
	{
		isTrim = true;
		trim = tr;
	}

	public boolean isTrim()
	{
		return isTrim;
	}

	public double getTrimValue()
	{
		return trim;
	}

	public void setTop(String top)
	{
		this.Top = Integer.parseInt(top);
	}
	
	public int getTop()
	{
		return Top;
	}

//	public boolean isFst()
//	{
//		return isFst;
//	}

	public int getMode()
	{
		return mode;
	}

	protected static int BETA = 0;
	protected static int FRQ = 1;

	private ArrayList<String> md;
	private boolean isGZ = false;
	private boolean isQT = true;
	private boolean isVerbose = false;
	private boolean isVerboseGZ = false;

	private boolean isClean = false;
	private boolean isRapid = false;
	private boolean isRandom = false;
	private boolean isTrim = false;
	private double trim = 0;

	private boolean isBeta = true;
	private boolean isFrq = false;
	private int mode = BETA;

	private double Me = 30000;
	private int Top = 0;

	private int chr = 0;
	private boolean chrFlag = false;

	private double[] qtSize;
	private double[] ccSize;
	private String[] field = {"snp", "chr", "bp", "beta", "or", "se", "p", "a1", "a2"};
	
}
