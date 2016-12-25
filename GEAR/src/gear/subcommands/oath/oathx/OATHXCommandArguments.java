package gear.subcommands.oath.oathx;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OATHXCommandArguments extends CommandArguments 
{


	public void setChr(String c)
	{
		this.chr = Integer.parseInt(c);
		if (this.chr < 1)
		{
			Logger.printUserLog("Chromosome should be greater than 0.\n GEAR quitted.");
			System.exit(1);
		}
		this.chrFlag = true;
	}

	public int getChr()
	{
		return this.chr;
	}

	public boolean isChrFlagOn()
	{
		return this.chrFlag;
	}

	public void setKeeFile(String kFile) 
	{
		FileUtil.exists(kFile);
		keepFile = kFile;
	}
	
	public String getKeepFile()
	{
		return keepFile;
	}
	
	public void setMAF(String mf) 
	{
		maf = Double.parseDouble(mf);
	}
	
	public double getMAF()
	{
		return maf;
	}

	public void setN(String ns)
	{	
		n = Integer.parseInt(ns);
		if (n < 1)
		{
			Logger.printUserError("--n should be greater than 1");
			System.exit(1);
		}
	}
	
	public int getN()
	{
		return n;
	}

	public void setNSSBatch(String batch) 
	{
		FileUtil.exists(batch);
		md = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "MetaBatch");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			FileUtil.exists(tokens[0]);
			md.add(tokens[0]);
		}
		batchFile = batch;
	}

	public String getBatchFile()
	{
		return batchFile;
	}

	public String[] getNSSFile() 
	{
		return md.toArray(new String[0]);
	}

	public void setCMFile(String cm) 
	{
		FileUtil.exists(cm);
		cmFile = cm;
	}
	
	public String getCMFile()
	{
		return cmFile;
	}

	public void setVerbose() 
	{
		verbose = true;
	}
	
	public boolean isVerbose()
	{
		return verbose;
	}

	private int chr;
	private int n;
	private boolean chrFlag = false;
	private String cmFile = null;
	private String batchFile = null;
	private ArrayList<String> md;

	private String keepFile;

	private double maf = 0.05;
	private boolean verbose = false;

}
