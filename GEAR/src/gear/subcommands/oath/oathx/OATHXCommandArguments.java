package gear.subcommands.oath.oathx;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OATHXCommandArguments extends CommandArguments 
{
	
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

	private int n;
	private String cmFile = null;
	private String batchFile = null;
	private ArrayList<String> md;

	private double maf = 0.05;
	private boolean verbose = false;

}
