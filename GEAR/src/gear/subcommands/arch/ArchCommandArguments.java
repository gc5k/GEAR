package gear.subcommands.arch;

import java.util.ArrayList;

import gear.gwassummary.GWASConstant;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class ArchCommandArguments extends CommandArguments
{
	public void setWindow(String k)
	{
		window = Double.parseDouble(k);
		windowFlag = true;
	}

	public double getWindow()
	{
		return window;
	}
	
	public boolean isWindow()
	{
		return windowFlag;
	}

	public void setLDFile(String f)
	{
		ldrFile = new String(f + ".gz");
		FileUtil.exists(ldrFile);
		ldsnpFile = new String(f + ".snp");
		FileUtil.exists(ldsnpFile);
	}
	
	public String getLDRFile()
	{
		return ldrFile;
	}
	
	public String getLDSNPFile()
	{
		return ldsnpFile;
	}

	public void setQtMetaBatch(String batch)
	{
		FileUtil.exists(batch);
		md = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "Qt Meta Batch");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			FileUtil.exists(tokens[0]);
			md.add(tokens[0]);
		}
		reader.close();
		isQT = true;
	}

	public void setCCMetaBatch(String batch)
	{
		FileUtil.exists(batch);
		md = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "CC Meta Batch");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			FileUtil.exists(tokens[0]);
			md.add(tokens[0]);
		}
		reader.close();
		isQT = false;
	}

	public String[] getMetaFile() 
	{
		return md.toArray(new String[0]);
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

	public boolean isQT()
	{
		return isQT;
	}

	public boolean isGZ()
	{
		return isGZ;
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

	public void setExtract(String extract)
	{
		this.extract = extract;
		FileUtil.exists(extract);
	}
	
	public String getExtractFile()
	{
		return extract;
	}

	private String extract;
	private boolean chrFlag = false;
	private int chr = 0;

	private boolean isGZ = false;
	private boolean isQT = true;
	private String ldrFile = null;
	private String ldsnpFile = null;

	private ArrayList<String> md;
	private String[] field = {"snp", "chr", "bp", "beta", "or", "se", "p", "a1", "a2"};

	private double window = 1;
	private boolean windowFlag = true;
}
