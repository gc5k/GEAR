package gear.util.structure;

import gear.ConstValues;

public class MAF
{
	private String chr;
	private String snp;
	private char allele1;
	private char allele2;
	private double maf;
	private double nchr;
	
	public MAF(String chr, String snp, char allele1, char allele2, double maf, double nchr)
	{
		this.chr = chr;
		this.snp = snp;
		this.allele1 = allele1;
		this.allele2 = allele2;
		this.maf = maf;
		this.nchr = nchr;
	}

	public String getChr()
	{
		return chr;
	}

	public String getSNP()
	{
		return snp;
	}

	public char getAllele1()
	{
		return allele1;
	}

	public char getAllele2()
	{
		return allele2;
	}

	public double getMAF()
	{
		return maf;
	}

	public double getNChr()
	{
		return nchr;
	}
	
	public static MAF next(gear.util.BufferedReader reader)
	{
		String[] tokens = reader.readTokens(6);
		
		if (tokens == null)
		{
			return null;
		}
		
		String chr = tokens[0];
		String snp = tokens[1];
		
		if (tokens[2].length() != 1)
		{
			reader.reportFormatError("'" + tokens[2] + "' is not a valid allele. Alleles must be single characters.");
		}
		char allele1 = tokens[2].charAt(0);
		
		if (tokens[3].length() != 1)
		{
			reader.reportFormatError("'" + tokens[3] + "' is not a valid allele. Alleles must be single characters.");
		}
		char allele2 = tokens[3].charAt(0);

		double maf = Double.NaN;
		try
		{
			maf = ConstValues.isNA(tokens[4]) ? Double.NaN : Double.parseDouble(tokens[4]);
		}
		catch (NumberFormatException e)
		{
			reader.reportFormatError("'" + tokens[4] + "' is not a valid MAF value. An MAF value must be a floating point number.");
		}
		
		if (maf != Double.NaN && (maf < 0.0f || 1.0f < maf))
		{
			reader.reportFormatError("'" + tokens[4] + "' is not a valid MAF value. MAF must between 0 and 1 inclusively.");
		}
		
		double nchr = Double.NaN;
		try
		{
			nchr = Double.parseDouble(tokens[5]);
		}
		catch (NumberFormatException e)
		{
			reader.reportFormatError("'" + tokens[5] + "' is not a valid NCHR value. An NCHR value must be a floating point number.");
		}
		
		return new MAF(chr, snp, allele1, allele2, maf, nchr);
	}
}
