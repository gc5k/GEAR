package gear.subcommands.profile;

public class ScaleMaf
{
	protected ScaleMaf (String snp, char a1, double maf)
	{
		this.snp = snp;
		this.a1 = a1;
		this.maf = maf;
	}

	protected String getSNP()
	{
		return snp;
	}
	
	protected char getA1()
	{
		return a1;
	}
	
	protected double getMaf()
	{
		return maf;
	}
	private String snp;
	private char a1;
	private double maf;
}
