package gear.subcommands.metapc;

import java.util.ArrayList;

import gear.subcommands.metapc.freader.FConstant;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class MetaPCCommandArguments extends CommandArguments
{
	public void setMetaBatch(String batch)
	{
		FileUtil.exists(batch);
		md = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "MetaBatch");

//		Logger.printUserLog("Checking the summary statistic files...");
		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			FileUtil.exists(tokens[0]);
			md.add(tokens[0]);
		}
		reader.close();
//		Logger.printUserLog("Found all of " + md.size() + " files.");
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
		reader.close();

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
		reader.close();

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
		field[FConstant.SNP] = k[0];
		field[FConstant.CHR] = k[1];

		field[FConstant.A1] = k[2];
		field[FConstant.A2] = k[3];
		field[FConstant.Fvalue] = k[4];
	

//		if(k.length > 5)
//		{
//			field[FConstant.BP] = k[5];
//			Logger.printUserLog("[INFO] The keyword for 'BP' (base pair) is set to " + field[FConstant.BP]);
//			keyLen=6;
//		}
//
//		if(k.length > 6)
//		{
//			field[FConstant.BP] = k[6];
//			Logger.printUserLog("[INFO] The keyword for 'N' (sample size) is set to " + field[FConstant.N]);
//			keyLen=7;
//		}
	}

	public String[] getKeys()
	{
		return field;
	}

//	public void setVerbose()
//	{
//		isVerbose = true;
//	}
//
//	public boolean isVerbose()
//	{
//		return isVerbose;
//	}

//	public double getMe()
//	{
//		return Me;
//	}
//
//	public void setNoWeight()
//	{
//		isNoWeight = true;
//	}
//	
//	public boolean isNoWeight()
//	{
//		return isNoWeight;
//	}

	public double getNe()
	{
		return Ne;
	}

//	public void setMe(String me)
//	{
//		Me = Double.parseDouble(me);
//	}

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

	public void setBeta()
	{
		isBeta = true;
	}
	
	public boolean isBeta()
	{
		return isBeta;
	}
	
	public void keepATGC()
	{
		isKeepATGC = true;
	}

	public boolean isKeepATGC()
	{
		return isKeepATGC;
	}

//	public void setRapid()
//	{
//		isRapid = true;
//	}
//
//	public boolean isRapid()
//	{
//		return isRapid;
//	}

//	public void setMeFrac(double mefrac)
//	{
//		meFrac = mefrac;
//	}
//
//	public double getMeFrac()
//	{
//		return meFrac;
//	}

//	public void setTop(String top)
//	{
//		this.Top = Integer.parseInt(top);
//	}
//
//	public int getTop()
//	{
//		return Top;
//	}

//	public boolean isFrq()
//	{
//		return isFrq;
//	}
//
//	public int getMode()
//	{
//		return mode;
//	}

	public String toString()
	{
		StringBuffer str = new StringBuffer();
		str.append("\nThe keyword for 'SNP' is set to " + field[FConstant.SNP] + "\n");
		str.append("The keyword for 'CHR' is set to " + field[FConstant.CHR] + "\n");
		str.append("The keyword for 'A1' (reference allele) is set to " + field[FConstant.A1] + "\n");
		str.append("The keyword for 'A2' (the other allele) is set to " + field[FConstant.A2] + "\n");
		str.append("The keyword for 'Value' (reference value for A1) is set to " + field[FConstant.Fvalue] + "\n");

//		if(keyLen > 5)
//		{
//			str.append("The keyword for 'BP' (base pair) is set to " + field[FConstant.BP] + "\n");
//		}
//
//		if(keyLen > 6)
//		{
//			str.append("The keyword for 'N' (sample size) is set to " + field[FConstant.N] + "\n");
//		}

		return str.toString();
	}

//	protected static int FRQ = 1;
//	private int mode = FRQ;

//	public boolean isFrq = true;

	private ArrayList<String> md;
	private boolean isGZ = false;
	private boolean isQT = true;
//	private boolean isVerbose = false;

//	private boolean isNoWeight = false;
	private double Ne = 500;
//	private boolean isRapid = false;
//	private int keyLen = 5;

//	private double Me = 30000;
//	private double meFrac = 0;
//	private int Top = 0;

	private int chr = 0;
	private boolean chrFlag = false;

	private boolean isBeta = false;
	private boolean isKeepATGC = false;

	private double[] qtSize;
	private double[] ccSize;
	private String[] field = {"snp", "chr", "a1", "a2", "value"};//, "or", "se", "p", };


}
